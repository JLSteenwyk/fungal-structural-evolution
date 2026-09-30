#!/usr/bin/env python3
"""Full synthetic native union: source direction, exclusions, errors and false exports."""
import csv
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
import numpy as np
from duplication_alignment_inputs import render_ca
from duplication_alignment_numeric_readback import load_pdb
from duplication_alignment_numeric_diagnostic import check_alignment
from assess_reference_alignment_geometry import geometry
from run_duplication_alignments import run_job
from union_reference_measurements import exclusions
from run_ortholog_pair_guide_comparison import sha


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2)+'\n')


def table(path, rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def main():
    assert exclusions('aligned',dict(rmsd_status='outside_printed_rounding',aligned_length='2'),dict(geometry_status='degenerate_at_numeric_tolerance')) == ['rmsd_discrepancy','fewer_than_three_pairs','nonunique_rotation']
    assert exclusions('timeout',None,None)==['timeout']
    with tempfile.TemporaryDirectory(prefix='reference-union-fixture-') as temp:
        root=Path(temp);seq='ACDEFGHIKLMNPQRSTVWY';keys=dict(A=('A',6),B=('B',10),C=('C',6),D=('D',10))
        inputs={}
        for name,end in keys.items():
            xyz=[[float(i),0.,0.] for i in range(len(seq))] if name=='C' else [[3*math.cos(i),3*math.sin(i),i] for i in range(len(seq))]
            for mask in ['full','plddt70']:
                rec=dict(status='validated',sequence=seq,ca_xyz=xyz,ca_plddt=[90.]*len(seq))
                if name=='C' and mask=='plddt70':rec['ca_plddt']=[90.,90.]+[0.]*(len(seq)-2)
                blob,letters,positions=render_ca(rec,threshold=70 if mask=='plddt70' else None)
                row=dict(sequence=letters,status='too_few_retained_residues',original_positions=positions,retained_residues=len(positions))
                if len(positions)>=3:
                    path=root/(name+'-'+mask+'.pdb');path.write_bytes(blob);row.update(status='ready',path=str(path),sha256=sha(path))
                inputs[(*end,mask)]=row
        pairs=[]
        for a,b in [('A','B'),('A','C'),('B','C'),('A','D')]:
            ends=[keys[a],keys[b]];pair=hashlib.sha256(json.dumps(sorted(ends),separators=(',',':')).encode()).hexdigest()
            pairs.append(dict(pair_key=pair,model_a=ends[0][0],version_a=ends[0][1],model_b=ends[1][0],version_b=ends[1][1],measurement_disposition='native_measurement_required' if a=='B' or b=='D' else 'pending_actual_input_result_and_numeric_reuse_checks'))
        native=root/'native';native.mkdir();diagnostic=root/'diagnostic';diagnostic.mkdir();geo=root/'geometry';geo.mkdir()
        npth=root/'native-plan.json';dpth=root/'diagnostic-plan.json';gpth=root/'geometry-plan.json';qpth=root/'quaternion-plan.json';qp=root/'quaternion.json'
        exe='/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign'
        npconf=dict(output=str(native),usalign=exe,options=['-mol','prot','-mm','0','-outfmt','0','-ter','2'],per_pair_timeout_seconds=10,pins={exe:sha(exe)})
        write(npth,npconf);write(dpth,dict(output=str(diagnostic),source_plan=str(npth)));write(gpth,dict(output=str(geo),diagnostic_plan=str(dpth)));write(qpth,dict(output=str(qp),source_plan=str(gpth)))
        checkpoints=[];numerics=[];geometries=[];reused=[]
        for index,p in enumerate(pairs):
            ends=[(p['model_a'],p['version_a']),(p['model_b'],p['version_b'])]
            source_ends=ends[::-1] if index==0 else ends
            oldfolder=root/('old'+str(index));oldfolder.mkdir();oldplan=dict(npconf,output=str(oldfolder));oldpath=root/('old-plan'+str(index)+'.json');write(oldpath,oldplan)
            config=oldplan if index<2 else npconf;ph=sha(oldpath if index<2 else npth)
            for mask in ['full','plddt70']:
                for source_order in [0,1]:
                    directed=source_ends if source_order==0 else source_ends[::-1]
                    path,status=run_job((p['pair_key'],*directed,mask,source_order),inputs,config,ph,'fixture-input-bundle')
                    raw=json.loads(path.read_text())
                    if index==3 and (mask=='plddt70' or source_order==1):
                        raw.pop('metrics',None)
                        if mask=='full':raw.update(status='parse_error',error='synthetic parse disposition',returncode=0)
                        elif source_order==0:raw.update(status='timeout',timeout_seconds=10)
                        else:raw.update(status='native_error',returncode=1)
                        write(path,raw);status=raw['status']
                    n=g=None
                    if status=='aligned':
                        coords=[load_pdb(inputs[(*end,mask)]) for end in directed]
                        n=dict(pair_key=p['pair_key'],mask=mask,order=str(source_order),**{k:str(v) for k,v in check_alignment(raw,*coords).items()})
                        strings=[raw['metrics'][f'alignment_{s}'] for s in ['left','right']];nongap=[np.array(list(s))!='-' for s in strings];paired=nongap[0]&nongap[1];ix=[(np.cumsum(a)-1)[paired] for a in nongap]
                        g=dict(pair_key=p['pair_key'],mask=mask,order=str(source_order),rmsd_status=n['rmsd_status'],**{k:str(v) for k,v in geometry(*[c[1][i] for c,i in zip(coords,ix)]).items()})
                        if index==1 and mask=='full' and source_order==0:
                            n['rmsd_status']=g['rmsd_status']='outside_printed_rounding'  # Explicit synthetic quarantine fixture.
                    order=1-source_order if index==0 else source_order
                    if index>=2:
                        checkpoints.append(dict(path=str(path.relative_to(native)),sha256=sha(path),status=status))
                        if n is not None:numerics.append(n);geometries.append(g)
                        r=dict(selected_source='',reuse_status='new_native_measurement_pending',source_checkpoint=None,source_checkpoint_sha256=None,source_order=None,source_native_status=None,source_numeric=None,source_geometry=None,numerical_exclusion_reasons=[],numerical_usable=False)
                    else:
                        why=exclusions(status,n,g)
                        r=dict(selected_source='primary_expanded_completed' if index==0 else 'reference_old',reuse_status='verified_identical_input_checkpoint_and_retained_disposition',source_checkpoint=str(path),source_checkpoint_sha256=sha(path),source_order=source_order,source_native_status=status,source_numeric=n,source_geometry=g,numerical_exclusion_reasons=why,numerical_usable=not why)
                    reused.append(dict(**{k:v for k,v in p.items() if k!='measurement_disposition'},mask=mask,order=order,**r))
        table(native/'full_reference_work_partition.tsv',pairs);table(diagnostic/'full_reference_work_partition.tsv',pairs);table(native/'checkpoint_manifest.tsv',checkpoints);table(diagnostic/'numeric_readback.tsv',numerics);table(geo/'alignment_geometry.tsv',geometries)
        counts=dict(Counter(Path(r['path']).stem.rsplit('-',2)[1]+':'+r['status'] for r in checkpoints));gcounts=dict(Counter(r['mask']+':'+r['geometry_status'] for r in geometries))
        write(native/'receipt.json',dict(status='complete_reference_alignment_dispositions_pending_readback',plan_sha256=sha(npth),counts=counts,directed_dispositions=8,distinct_model_pairs=2,full_reference_pairs=4,existing_catalog_pairs_pending_reuse=2,upstream_bindings={},input_bundle_sha256='fixture-input-bundle',artifacts={n:sha(native/n) for n in ['full_reference_work_partition.tsv','checkpoint_manifest.tsv']}))
        write(diagnostic/'receipt.json',dict(status='complete_reference_alignment_rmsd_diagnostic_not_scientific_acceptance',plan_sha256=sha(dpth),producer_receipt_sha256=sha(native/'receipt.json'),numerically_checked_alignments=len(numerics),directed_dispositions=8,existing_catalog_pairs_pending_reuse=2,counts=counts,rmsd_status_counts=dict(Counter(r['rmsd_status'] for r in numerics)),artifacts={n:sha(diagnostic/n) for n in ['full_reference_work_partition.tsv','numeric_readback.tsv']}))
        write(geo/'receipt.json',dict(status='complete_reference_alignment_geometry_pending_independent_readback',plan_sha256=sha(gpth),diagnostic_receipt_sha256=sha(diagnostic/'receipt.json'),alignments=len(geometries),counts=gcounts,artifacts={'alignment_geometry.tsv':sha(geo/'alignment_geometry.tsv')}))
        write(qp,dict(status='passed_full_reference_geometry_readback',plan_sha256=sha(qpth),producer_receipt_sha256=sha(geo/'receipt.json'),alignments_checked=len(geometries),counts=gcounts))
        na=root/'native-archive.json';write(na,dict(services=[{'synthetic_stub':True}]*4,summary=dict(directed_dispositions=8),source_hashes={str(native/'receipt.json'):sha(native/'receipt.json'),str(qp):sha(qp)}))
        nc=root/'native-closure.json';write(nc,dict(status='complete_verified_new_reference_native_alignment_numeric_geometry',new_pairs=2,full_pairs=4,full_hash_archive=str(na),full_hash_archive_sha256=sha(na)))
        rr=root/'reuse';rr.mkdir();(rr/'reference_reuse_dispositions.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in sorted(reused,key=lambda r:(r['pair_key'],r['mask'],r['order']))))
        rcounts=dict(new_native_measurement_pending=8,verified_identical_input_checkpoint_and_retained_disposition=8)
        write(rr/'receipt.json',dict(status='complete_full_reference_reuse_qualification_pending_independent_readback',counts=rcounts,artifacts={'reference_reuse_dispositions.jsonl':sha(rr/'reference_reuse_dispositions.jsonl')}))
        ra=root/'reuse-archive.json';write(ra,dict(summary=dict(directed_dispositions=16),services=[{'synthetic_stub':True}]*2))
        rp=root/'reuse-proof.json';write(rp,dict(status='passed_full_reference_input_checkpoint_numeric_reuse_readback',producer_receipt_sha256=sha(rr/'receipt.json'),counts=rcounts))
        rs=root/'reuse-plan.json';write(rs,dict(target_alignment_plan=str(npth)))
        rc=root/'reuse-closure.json';write(rc,dict(status='complete_verified_full_reference_input_checkpoint_numeric_reuse',full_reference_pairs=4,directed_dispositions=16,counts=rcounts,full_hash_archive=str(ra),full_hash_archive_sha256=sha(ra),producer_receipt=str(rr/'receipt.json'),producer_receipt_sha256=sha(rr/'receipt.json'),independent_readback=str(rp),independent_readback_sha256=sha(rp),source_plan=str(rs),source_plan_sha256=sha(rs)))
        plan=root/'plan.json';out=root/'union';write(plan,dict(full_pairs=4,new_pairs=2,native_plan=str(npth),native_diagnostic_plan=str(dpth),native_geometry_plan=str(gpth),native_geometry_readback_plan=str(qpth),native_readback=str(qp),native_completion=str(nc),reuse_completion=str(rc),output=str(out),pins={}))
        commands=[[sys.executable,'scripts/union_reference_measurements.py','--plan',str(plan)], [sys.executable,'scripts/readback_reference_measurement_union.py','--plan',str(plan),'--output',str(root/'proof.json')]]
        for command in commands:
            run=subprocess.run(command,capture_output=True,text=True)
            if run.returncode:raise RuntimeError(run.stderr+run.stdout)
        receipt=json.loads((out/'receipt.json').read_text());assert receipt['directed_dispositions']==16
        original=(out/'reference_measurement_dispositions.jsonl').read_text();rejected=[]
        for label in ['cleared_exclusion','wrong_source_order','reversed_directed_endpoints','changed_native_metric','changed_checkpoint_hash','missing_unavailable_state']:
            rows=[json.loads(l) for l in original.splitlines()]
            if label=='cleared_exclusion':
                r=next(r for r in rows if r['source_native_status']=='aligned' and not r['numerical_usable']);r['numerical_usable']=True;r['numerical_exclusion_reasons']=[]
            elif label=='wrong_source_order':rows[0]['source_order']=1-rows[0]['source_order']
            elif label=='reversed_directed_endpoints':rows[0]['directed_endpoints'].reverse()
            elif label=='changed_native_metric':next(r for r in rows if r['native_metrics'])['native_metrics']['rmsd']=999
            elif label=='changed_checkpoint_hash':rows[0]['source_checkpoint_sha256']='0'*64
            else:rows.pop(next(i for i,r in enumerate(rows) if r['source_native_status']=='input_unavailable'))
            (out/'reference_measurement_dispositions.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows));receipt['artifacts']['reference_measurement_dispositions.jsonl']=sha(out/'reference_measurement_dispositions.jsonl');write(out/'receipt.json',receipt)
            run=subprocess.run(commands[1][:-1]+[str(root/(label+'.json'))],capture_output=True,text=True);assert run.returncode!=0,label;rejected.append(label)
        print(json.dumps(dict(status='passed_full_synthetic_native_reference_measurement_union_checks',full_pairs=4,directed_states=16,rejected_rehashed_exports=rejected,scope='Actual native USalign outputs with versions6/10, reversed original primary directions, both masks, unavailable inputs, degenerate rotation and explicit synthetic RMSD/error/parse/timeout dispositions. Original qualification/numeric/quaternion/journal proof stubs are fixture data only; not production qualification or a pilot. All three exclusion causes tested directly.'),indent=2))


if __name__=='__main__':main()
