#!/usr/bin/env python3
"""Check all reference order correspondences and quantiles independently."""
import csv,json,hashlib,math
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    root=Path('results/structural_comparisons/duplication-reference-order-sensitivity-20260927-v1');r=json.loads((root/'receipt.json').read_text())
    assert r['status']=='complete_full_reference_order_sensitivity_pending_independent_readback'
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    source=Path('results/structural_comparisons/duplication-reference-usable-orders-20260927-v1');sr=json.loads((source/'receipt.json').read_text())
    assert sha(source/'receipt.json')==r['source_receipt_sha256'] and sha(source/'pair_mask_order_summary.tsv')==sr['artifacts']['pair_mask_order_summary.tsv']
    original={(x['pair_key'],x['mask']):x for x in rows(source/'pair_mask_order_summary.tsv') if x['order_summary_status']=='both_orders_numerically_usable'}
    native=Path('results/structural_comparisons/duplication-reference-alignments-20260926-v1');nr=json.loads((native/'receipt.json').read_text())
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
    result=dict(status='passed_full_reference_order_sensitivity_readback',pair_mask_rows=len(seen),native_mappings_checked=2*len(seen),quantile_rows=len(checked),mapping_counts=r['mapping_counts'],source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(Path(__file__)),scope='All correspondence sets independently reconstructed with cumulative-index arrays; all intersection/union/status/metric fields and common cohorts checked; quantiles verified by sorted interpolation. No alignment rerun or biological qualification.')
    with Path('metadata/duplication_reference_order_sensitivity_readback_20260927.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
