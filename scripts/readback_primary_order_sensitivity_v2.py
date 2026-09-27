#!/usr/bin/env python3
"""Check all primary order correspondences and quantiles independently."""
import argparse,csv,json,math,subprocess,time
import psutil
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np


from screen_duplication_domain_alignment_coverage import sha


def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for path,digest in plan['pins'].items():assert sha(path)==digest,path
    verify();dep=plan['producer']
    while True:
        try:
            proc=psutil.Process(dep['pid'])
            if proc.create_time()!=dep['created'] or proc.status()==psutil.STATUS_ZOMBIE:break
            assert proc.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        print('waiting_for_verified_prerequisite',dep['pid'],flush=True);time.sleep(30)
    terminal=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert terminal==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),terminal
    verify()
    root=Path(plan['source']);r=json.loads((root/'receipt.json').read_text())
    assert r['plan_sha256']==sha(plan['source_plan'])
    assert r['source_readback_sha256']==sha(plan['orders_audit'])
    proof=json.loads(Path(plan['orders_audit']).read_text())
    assert proof['status']=='passed_full_primary_usable_order_summary_readback'
    assert proof['source_receipt_sha256']==r['source_receipt_sha256']
    assert r['status']=='complete_full_primary_order_sensitivity_pending_independent_readback'
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    source=Path(plan['orders']);sr=json.loads((source/'receipt.json').read_text())
    assert sha(source/'receipt.json')==r['source_receipt_sha256'] and sha(source/'pair_mask_order_summary.tsv')==sr['artifacts']['pair_mask_order_summary.tsv']
    original={(x['pair_key'],x['mask']):x for x in rows(source/'pair_mask_order_summary.tsv') if x['order_summary_status']=='both_orders_numerically_usable'}
    native=Path(plan['native']);nr=json.loads((native/'receipt.json').read_text())
    assert sha(native/'receipt.json')==r['native_receipt_sha256']
    assert sha(native/'checkpoint_manifest.tsv')==nr['artifacts']['checkpoint_manifest.tsv']
    hashes={x['path']:x['sha256'] for x in rows(native/'checkpoint_manifest.tsv')}
    common={p for p,m in original if m=='full'} & {p for p,m in original if m=='plddt70'}
    seen=set();groups=defaultdict(list);statuses=defaultdict(Counter)
    data=rows(root/'pair_mask_sensitivity.tsv')
    metrics=[k for k in data[0] if k.endswith('_absolute_difference')]
    for row in data:
        pair,mask=row['pair_key'],row['mask'];key=pair,mask
        assert key not in seen and key in original;seen.add(key);sets=[]
        for order in [0,1]:
            rel=f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json';path=native/rel;assert sha(path)==hashes[rel]
            record=json.loads(path.read_text());strings=[record['metrics'][name] for name in ['alignment_left','alignment_right']]
            masks=[np.array(list(x))!='-' for x in strings];shared=masks[0]&masks[1];indices=[(np.cumsum(x)-1)[shared] for x in masks]
            if order==1:indices.reverse()
            sets.append(set(zip(*indices)))
        a,b=sets;status='identical' if a==b else 'different_same_count' if len(a)==len(b) else 'different_count'
        assert row['mapping_status']==status
        for field,value in [('order0_pairs',len(a)),('order1_pairs',len(b)),('shared_pairs',len(a&b)),('union_pairs',len(a|b))]:assert int(row[field])==value
        assert float(row['mapping_jaccard'])==len(a&b)/len(a|b)
        assert row['in_both_mask_cohort']==str(pair in common)
        cohorts=['full_all'] if mask=='full' else []
        if pair in common:cohorts.append('full_common' if mask=='full' else 'plddt70_common')
        for metric in metrics:
            name=metric.removesuffix('_absolute_difference');expected=abs(float(original[key]['order0_'+name])-float(original[key]['order1_'+name]));assert float(row[metric])==expected
        for cohort in cohorts:
            statuses[cohort][status]+=1
            for metric in metrics+['mapping_jaccard']:groups[cohort,metric].append(float(row[metric]))
    assert seen==set(original) and len(seen)==r['pair_mask_rows'] and len(common)==r['common_mask_pairs']
    assert {k:dict(v) for k,v in statuses.items()}==r['mapping_counts']
    checked=set()
    for row in rows(root/'quantiles.tsv'):
        key=row['cohort'],row['metric'];assert key not in checked;checked.add(key)
        values=sorted(groups[key]);assert len(values)==int(row['n'])
        for name,q in [('minimum',0),('median',.5),('q90',.9),('q95',.95),('q99',.99),('maximum',1)]:
            position=(len(values)-1)*q;i=math.floor(position);fraction=position-i
            expected=values[i]*(1-fraction)+values[min(i+1,len(values)-1)]*fraction
            assert math.isclose(float(row[name]),expected,rel_tol=1e-12,abs_tol=1e-12)
    assert checked==set(groups)
    verify()
    result=dict(plan_sha256=ph,prerequisite_terminal_state=terminal,status='passed_full_primary_order_sensitivity_readback',pair_mask_rows=len(seen),native_mappings_checked=2*len(seen),quantile_rows=len(checked),mapping_counts=r['mapping_counts'],source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(Path(__file__)),scope='All correspondence sets independently reconstructed with cumulative-index arrays; all intersection/union/status/metric fields and common cohorts checked; quantiles verified by sorted interpolation. No alignment rerun or biological qualification.')
    with Path(plan['output']).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
