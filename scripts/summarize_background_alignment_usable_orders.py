#!/usr/bin/env python3
"""Preserve every background_alignment pair and input order with explicit numerical exclusions."""
import argparse,csv,json,time,subprocess
from collections import Counter
from pathlib import Path
import psutil
from summarize_duplication_alignment_orders import summarize_pair,sha


def reasons(numeric,geometry):
    result=[]
    if numeric['rmsd_status']!='within_printed_rounding':result.append('rmsd_discrepancy')
    if int(numeric['aligned_length'])<3:result.append('fewer_than_three_pairs')
    if geometry['geometry_status']!='unique_at_numeric_tolerance':result.append('nonunique_rotation')
    return ';'.join(result)


def pair_summary(pair,mask,statuses,numeric,geometry):
    adjusted={};retained={};extra={}
    for order in [0,1]:
        raw=statuses[order];why=''
        if raw=='aligned':
            why=reasons(numeric[order],geometry[order])
            if not why:retained[order]=numeric[order]
        adjusted[order]='excluded_numerically' if why else raw
        extra[f'order{order}_native_status']=raw
        extra[f'order{order}_numerical_exclusion_reasons']=why
        extra[f'order{order}_geometry_status']=geometry[order]['geometry_status'] if raw=='aligned' else ''
        extra[f'order{order}_rmsd_status']=numeric[order]['rmsd_status'] if raw=='aligned' else ''
    row=summarize_pair(pair,mask,adjusted,retained)
    row['order_summary_status']=row['order_summary_status'].replace('aligned','numerically_usable')
    return dict(row,**extra)


def table(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args()
    plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();dep=plan['producer']
    while True:
        try:
            p=psutil.Process(dep['pid'])
            if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        print('waiting_for_full_background_alignment_geometry_verification',dep['pid'],flush=True);time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    if state!=dict(ActiveState='inactive',Result='success',ExecMainStatus='0'):raise ValueError('Geometry audit did not finish successfully: '+str(state))
    verify();root=Path(plan['native']);diagnostic=Path(plan['diagnostic']);geometry=Path(plan['geometry'])
    audit=json.loads(Path(plan['geometry_audit']).read_text());gr=json.loads((geometry/'receipt.json').read_text())
    assert audit['status']=='passed_full_background_geometry_readback' and audit['alignments_checked']==282120
    assert audit['producer_receipt_sha256']==sha(geometry/'receipt.json')
    assert audit['plan_sha256']==sha(plan['geometry_audit_plan'])
    dr=json.loads((diagnostic/'receipt.json').read_text());nr=json.loads((root/'receipt.json').read_text())
    assert dr['status']=='complete_background_alignment_rmsd_diagnostic_not_scientific_acceptance'
    assert gr['diagnostic_receipt_sha256']==sha(diagnostic/'receipt.json') and dr['producer_receipt_sha256']==sha(root/'receipt.json')
    for folder,receipt in [(root,nr),(diagnostic,dr),(geometry,gr)]:
        for name,h in receipt['artifacts'].items():assert sha(folder/name)==h
    statuses={}
    for row in table(root/'checkpoint_manifest.tsv'):
        pair,mask,order=Path(row['path']).stem.split('-');key=pair,mask,int(order)
        assert key not in statuses;statuses[key]=row['status']
    pairs={k[0] for k in statuses}
    assert len(pairs)==nr['distinct_model_pairs']==71450
    assert set(statuses)=={(p,m,o) for p in pairs for m in ['full','plddt70'] for o in [0,1]}
    numeric={};geos={}
    for target,path in [(numeric,diagnostic/'numeric_readback.tsv'),(geos,geometry/'alignment_geometry.tsv')]:
        for row in table(path):
            key=row['pair_key'],row['mask'],int(row['order']);assert key not in target;target[key]=row
    expected={k for k,v in statuses.items() if v=='aligned'}
    assert set(numeric)==set(geos)==expected and len(expected)==282120
    for key in expected:
        assert numeric[key]['aligned_length']==geos[key]['aligned_length'] and numeric[key]['rmsd_status']==geos[key]['rmsd_status']
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);counts=Counter();excluded=Counter()
    with (out/'pair_mask_order_summary.tsv').open('w') as f:
        writer=None
        for pair in sorted(pairs):
            for mask in ['full','plddt70']:
                row=pair_summary(pair,mask,{o:statuses[pair,mask,o] for o in [0,1]},{o:numeric[pair,mask,o] for o in [0,1] if (pair,mask,o) in numeric},{o:geos[pair,mask,o] for o in [0,1] if (pair,mask,o) in geos})
                if writer is None:writer=csv.DictWriter(f,list(row),delimiter='\t',lineterminator='\n');writer.writeheader()
                writer.writerow(row);counts[mask+':'+row['order_summary_status']]+=1
                for order in [0,1]:
                    why=row[f'order{order}_numerical_exclusion_reasons']
                    if why:excluded[why]+=1
    verify()
    result=dict(status='complete_background_alignment_numerically_usable_order_summary_pending_independent_readback',plan_sha256=ph,pairs=len(pairs),pair_mask_rows=2*len(pairs),counts=dict(counts),exclusion_reason_combinations=dict(excluded),geometry_audit_sha256=sha(plan['geometry_audit']),native_receipt_sha256=sha(root/'receipt.json'),diagnostic_receipt_sha256=sha(diagnostic/'receipt.json'),geometry_receipt_sha256=sha(geometry/'receipt.json'),artifacts={'pair_mask_order_summary.tsv':sha(out/'pair_mask_order_summary.tsv')},scope='Both input orders and every pair/mask retained. Numerically usable requires unchanged RMSD tolerance, at least three paired residues and unique rotation under the audited numerical criterion. Excluded metrics blank; reasons and native status retained. This is not coverage/confidence qualification, prediction accuracy, pair significance or scientific acceptance.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
    from readback_background_alignment_usable_orders import main as readback
    readback()
