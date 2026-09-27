"""Check every serialized primary geometry row against its native coordinates."""
import argparse
from collections import Counter
import csv
from functools import lru_cache
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import psutil
from duplication_alignment_numeric_readback import load_pdb
from readback_reference_alignment_geometry_v2 import verify_row
from screen_duplication_domain_alignment_coverage import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for p,h in plan['pins'].items(): assert sha(p)==h,p
    verify();launch=json.loads(Path(plan['producer_launch']).read_text())
    while True:
        try:
            process=psutil.Process(launch['pid'])
            if process.create_time()!=launch['created'] or process.status()==psutil.STATUS_ZOMBIE: break
            assert process.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess: break
        time.sleep(30)
    state=subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True)
    assert dict(l.split('=',1) for l in state.splitlines())==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    verify();gp=json.loads(Path(plan['producer_plan']).read_text());root=Path(gp['output'])
    assert sha(plan['producer_plan'])==launch['plan_sha256']
    rp=root/'receipt.json';rh=sha(rp);receipt=json.loads(rp.read_text())
    assert receipt['status']=='complete_primary_alignment_geometry_with_quaternion_checks'
    assert receipt['plan_sha256']==sha(plan['producer_plan'])
    for p,h in receipt['source_sha256'].items(): assert sha(p)==h,p
    table=root/'alignment_geometry.tsv';assert sha(table)==receipt['artifacts'][table.name]
    al=json.loads(Path(gp['audit_launch']).read_text());ap=json.loads(Path(al['plan']).read_text())
    assert sha(al['plan'])==al['plan_sha256']
    sp=json.loads(Path(ap['source_plan']).read_text());native=Path(sp['output'])
    audited=Path(ap['output']);ar=json.loads((audited/'receipt.json').read_text())
    assert ar['status']=='passed_full_duplication_alignment_mapping_rmsd_identity_readback'
    assert ar['producer_receipt_sha256']==sha(native/'receipt.json')
    manifest={}
    for row in csv.DictReader((native/'checkpoint_manifest.tsv').open(),delimiter='\t'):
        assert row['path'] not in manifest
        manifest[row['path']]=row
    expected={p for p,r in manifest.items() if r['status']=='aligned'}
    assert len(manifest)==receipt['directed_dispositions']
    assert len(expected)==receipt['alignments']==ar['numerically_checked_alignments']
    inputs={}
    for line in (Path(sp['inputs'])/'inputs.jsonl').open():
        r=json.loads(line);key=r['model_id'],r['version'],r['mask']
        assert key not in inputs;inputs[key]=r
    @lru_cache(maxsize=128)
    def coordinates(key): return load_pdb(inputs[key])
    seen=set();counts=Counter();near=0;maximum=0.
    for row in csv.DictReader(table.open(),delimiter='\t'):
        pair,mask,order=row['pair_key'],row['mask'],int(row['order'])
        rel=f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json'
        assert rel in expected and rel not in seen;seen.add(rel)
        assert sha(native/rel)==manifest[rel]['sha256']
        record=json.loads((native/rel).read_text())
        assert (record['pair_key'],record['mask'],record['order'],record['status'])==(pair,mask,order,'aligned')
        coords=[coordinates((v['model_id'],v['version'],mask)) for v in record['inputs']]
        strings=[record['metrics']['alignment_left'],record['metrics']['alignment_right']]
        assert len(strings[0])==len(strings[1])
        nongap=[np.array(list(s))!='-' for s in strings];paired=nongap[0]&nongap[1]
        for s,c in zip(strings,coords): assert s.replace('-','')==c[0]
        indexes=[(np.cumsum(p)-1)[paired] for p in nongap]
        error,close=verify_row(row,*[c[1][i] for c,i in zip(coords,indexes)])
        assert float(row['quaternion_scaled_curvature_error'])==error
        assert row['near_zero_quaternion_gap']==str(close)
        counts[mask+':'+row['geometry_status']]+=1;near+=int(close);maximum=max(maximum,error)
        if len(seen)%10000==0: print('Serialized primary geometry checked',len(seen),'/',len(expected),flush=True)
    assert seen==expected and dict(counts)==receipt['counts']
    assert near==receipt['near_zero_quaternion_gaps']
    assert maximum==receipt['maximum_scaled_quaternion_curvature_error']
    verify();assert sha(rp)==rh and sha(table)==receipt['artifacts'][table.name]
    for p,h in receipt['source_sha256'].items(): assert sha(p)==h,p
    result=dict(status='passed_full_primary_geometry_serialized_readback',plan_sha256=ph,
        producer_receipt_sha256=rh,alignments_checked=len(seen),counts=dict(counts),
        near_zero_quaternion_gaps=near,maximum_scaled_quaternion_curvature_error=maximum,
        scope='Every serialized geometry field checked against hashed native PDB coordinates and alignment maps using alternate SVD and quaternion eigensystem. Complete successful alignment set checked; skipped dispositions remain upstream. Numeric identifiability does not establish prediction accuracy or biological asymmetry.')
    out=Path(plan['output']);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f: f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__': main()
