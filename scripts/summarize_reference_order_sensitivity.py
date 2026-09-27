#!/usr/bin/env python3
"""Quantify complete reference input-order sensitivity and residue-map changes."""
import csv,hashlib,json,math
from collections import Counter
from pathlib import Path
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def mapping(record):
    a=record['metrics']['alignment_left'];b=record['metrics']['alignment_right'];assert len(a)==len(b)
    i=j=0;pairs=set()
    for x,y in zip(a,b):
        if x!='-' and y!='-':pairs.add((i,j) if record['order']==0 else (j,i))
        i+=x!='-';j+=y!='-'
    assert len(pairs)==record['metrics']['aligned_length']
    return pairs


def main():
    root=Path('results/structural_comparisons/duplication-reference-usable-orders-20260927-v1')
    native=Path('results/structural_comparisons/duplication-reference-alignments-20260926-v1')
    proof=Path('metadata/duplication_reference_usable_orders_readback_20260927.json')
    p=json.loads(proof.read_text());r=json.loads((root/'receipt.json').read_text())
    assert p['status']=='passed_full_reference_usable_order_summary_readback' and p['source_receipt_sha256']==sha(root/'receipt.json')
    assert r['artifacts']['pair_mask_order_summary.tsv']==sha(root/'pair_mask_order_summary.tsv')
    nr=json.loads((native/'receipt.json').read_text());assert sha(native/'receipt.json')==r['native_receipt_sha256']
    manifest=native/'checkpoint_manifest.tsv';assert sha(manifest)==nr['artifacts'][manifest.name]
    hashes={x['path']:x['sha256'] for x in rows(manifest)}
    source=rows(root/'pair_mask_order_summary.tsv');selected=[x for x in source if x['order_summary_status']=='both_orders_numerically_usable']
    common={x['pair_key'] for x in selected if x['mask']=='full'} & {x['pair_key'] for x in selected if x['mask']=='plddt70'}
    metrics=['aligned_length','rmsd_recomputed','sequence_identity_exact','tm_a','tm_b','coverage_a','coverage_b','joint_plddt70_fraction']
    output=[]
    for ix,row in enumerate(selected,1):
        pair,mask=row['pair_key'],row['mask'];maps=[]
        for order in [0,1]:
            rel=f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json';path=native/rel;assert sha(path)==hashes[rel]
            record=json.loads(path.read_text());assert (record['pair_key'],record['mask'],record['order'],record['status'])==(pair,mask,order,'aligned')
            maps.append(mapping(record))
        shared=len(maps[0]&maps[1]);union=len(maps[0]|maps[1]);identical=maps[0]==maps[1]
        value=dict(pair_key=pair,mask=mask,in_both_mask_cohort=pair in common,order0_pairs=len(maps[0]),order1_pairs=len(maps[1]),shared_pairs=shared,union_pairs=union,mapping_jaccard=shared/union,mapping_status='identical' if identical else 'different_same_count' if len(maps[0])==len(maps[1]) else 'different_count')
        for metric in metrics:
            x,y=float(row['order0_'+metric]),float(row['order1_'+metric]);difference=abs(x-y)
            assert math.isfinite(difference) and difference==float(row[metric+'_order_absolute_difference'])
            value[metric+'_absolute_difference']=difference
        if identical:assert value['rmsd_recomputed_absolute_difference']<1e-10
        output.append(value)
        if ix%10000==0:print('Compared mappings',ix,'/',len(selected),flush=True)
    assert len(output)==62887 and len(common)==30346
    quantiles=[];counts={}
    for label,subset in [('full_all',[x for x in output if x['mask']=='full']),('full_common',[x for x in output if x['mask']=='full' and x['in_both_mask_cohort']]),('plddt70_common',[x for x in output if x['mask']=='plddt70' and x['in_both_mask_cohort']])]:
        counts[label]=dict(Counter(x['mapping_status'] for x in subset))
        for metric in [x+'_absolute_difference' for x in metrics]+['mapping_jaccard']:
            values=np.array([x[metric] for x in subset]);qs=np.quantile(values,[0,.5,.9,.95,.99,1],method='linear')
            quantiles.append(dict(cohort=label,metric=metric,n=len(values),minimum=qs[0],median=qs[1],q90=qs[2],q95=qs[3],q99=qs[4],maximum=qs[5]))
    out=Path('results/structural_comparisons/duplication-reference-order-sensitivity-20260927-v1');out.mkdir(parents=True,exist_ok=False)
    for name,data in [('pair_mask_sensitivity.tsv',output),('quantiles.tsv',quantiles)]:
        with (out/name).open('w') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    result=dict(status='complete_full_reference_order_sensitivity_pending_independent_readback',pair_mask_rows=len(output),common_mask_pairs=len(common),mapping_counts=counts,source_receipt_sha256=sha(root/'receipt.json'),source_readback_sha256=sha(proof),native_receipt_sha256=sha(native/'receipt.json'),script_sha256=sha(Path(__file__)),artifacts={x.name:sha(x) for x in out.iterdir()},scope='All numerically usable two-order comparisons; residue correspondences checked from hashed checkpoints and endpoint-normalized metric differences retained. Full-all and identical common pair cohorts separated. Differences can reflect different aligned residues, not uncertainty bounds or biological changes. No favorable order selected; no coverage/confidence qualification.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
