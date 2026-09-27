#!/usr/bin/env python3
"""Measure coordinate rank and rigid-fit curvature for every recorded whole-protein primary alignment after full diagnostic verification."""
import argparse,csv,json,time,subprocess
import psutil
from collections import Counter
from functools import lru_cache
from pathlib import Path
import numpy as np
from duplication_alignment_numeric_readback import load_pdb
from duplication_alignment_numeric_diagnostic import check_alignment
from screen_duplication_domain_alignment_coverage import sha


from assess_reference_alignment_geometry import geometry
from readback_reference_alignment_geometry_v2 import verify_row


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed plan')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Changed pinned input: '+p)
    verify()
    completion=json.loads(Path(plan['completion']).read_text())
    assert completion['status']=='complete_verified_primary_rmsd_diagnostic_and_short_geometry'
    for unit in completion['terminal_states']:
        state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),state
    for path,digest in completion['source_hashes'].items():assert sha(path)==digest
    identity=plan['producer']
    while True:
        try:
            process=psutil.Process(identity['pid'])
            if process.create_time()!=identity['created'] or process.status()==psutil.STATUS_ZOMBIE:break
            assert process.cmdline()==identity['cmdline']
        except psutil.NoSuchProcess:break
        print('waiting_for_complete_reference_diagnostic',identity['pid'],flush=True);time.sleep(30)
    verify();diagnostic=Path(plan['diagnostic']);dr=json.loads((diagnostic/'receipt.json').read_text())
    if dr['status']!='complete_primary_alignment_rmsd_diagnostic_not_scientific_acceptance':raise ValueError('Wrong diagnostic status')
    if sha(diagnostic/'numeric_readback.tsv')!=dr['artifacts']['numeric_readback.tsv']:raise ValueError('Changed diagnostic table')
    dp=json.loads(Path(plan['diagnostic_plan']).read_text());source=json.loads(Path(dp['source_plan']).read_text());root=Path(source['output'])
    if dr['plan_sha256']!=sha(plan['diagnostic_plan']) or dp['output']!=str(diagnostic):raise ValueError('Diagnostic provenance mismatch')
    if sha(root/'receipt.json')!=dr['producer_receipt_sha256']:raise ValueError('Changed alignment producer')
    producer=json.loads((root/'receipt.json').read_text());manifest=root/'checkpoint_manifest.tsv'
    if sha(manifest)!=producer['artifacts']['checkpoint_manifest.tsv']:raise ValueError('Changed checkpoint manifest')
    proofs={r['path']:r['sha256'] for r in csv.DictReader(manifest.open(),delimiter='\t')}
    inputs={}
    for spec in [dict(inputs=source['inputs'])]:
        folder=Path(spec['inputs']);ir=json.loads((folder/'receipt.json').read_text())
        if sha(folder/'inputs.jsonl')!=ir['artifacts']['inputs.jsonl']:raise ValueError('Changed PDB manifest')
        for line in (folder/'inputs.jsonl').open():
            row=json.loads(line);key=row['model_id'],row['version'],row['mask']
            if key in inputs:raise ValueError('Duplicate input identity')
            inputs[key]=row
    @lru_cache(maxsize=256)
    def coords(key):return load_pdb(inputs[key])
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);counts=Counter();seen=set();maxerr=0.
    fields=['pair_key','mask','order','rmsd_status']+list(geometry(np.zeros((3,3)),np.zeros((3,3))))
    with (out/'alignment_geometry.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for row in csv.DictReader((diagnostic/'numeric_readback.tsv').open(),delimiter='\t'):
            pair,mask,order=row['pair_key'],row['mask'],int(row['order']);key=pair,mask,order
            if key in seen:raise ValueError('Repeated alignment');
            seen.add(key);rel=f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json';path=root/rel
            if sha(path)!=proofs[rel]:raise ValueError('Changed checkpoint')
            record=json.loads(path.read_text())
            if record['status']!='aligned' or (record['pair_key'],record['mask'],record['order'])!=key:raise ValueError('Wrong checkpoint')
            left,right=[coords((r['model_id'],r['version'],mask)) for r in record['inputs']]
            numeric=check_alignment(record,left,right)
            for name,value in numeric.items():
                if name=='rmsd_status':
                    if row[name]!=value:raise ValueError('Changed diagnostic classification')
                elif float(row[name])!=float(value):raise ValueError('Changed diagnostic metric')
            x,y=record['metrics']['alignment_left'],record['metrics']['alignment_right'];i=j=0;ix=[];iy=[]
            for a,b in zip(x,y):
                if a!='-' and b!='-':ix.append(i);iy.append(j)
                i+=a!='-';j+=b!='-'
            g=geometry(left[1][ix],right[1][iy])
            verify_row(g,left[1][ix],right[1][iy])
            if g['aligned_length']!=int(row['aligned_length']):raise ValueError('Mapping length mismatch')
            w.writerow(dict(pair_key=pair,mask=mask,order=order,rmsd_status=row['rmsd_status'],**g))
            counts[mask+':'+g['geometry_status']]+=1;maxerr=max(maxerr,numeric['rmsd_rounding_error'])
            if len(seen)%10000==0:print('Geometry checked',len(seen),'/',dr['numerically_checked_alignments'],flush=True)
    if len(seen)!=dr['numerically_checked_alignments']:raise ValueError('Incomplete geometry cohort')
    verify()
    result=dict(status='complete_primary_diagnostic_geometry_pending_independent_readback',plan_sha256=ph,diagnostic_receipt_sha256=sha(diagnostic/'receipt.json'),alignments=len(seen),counts=dict(counts),maximum_rmsd_discrepancy=maxerr,artifacts={'alignment_geometry.tsv':sha(out/'alignment_geometry.tsv')},scope='All successful directed whole-protein primary alignment after full diagnostic verifications; native mappings and numeric diagnostic rechecked from hashed PDBs, coordinate singular spectra and optimal proper-rotation curvature recorded; every in-memory geometry row checked by the independent quaternion verifier. Serialized readback remains required. Numerical identifiability, not prediction accuracy or stability to coordinate/model uncertainty. Existing RMSD quarantine unchanged; no biological inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
