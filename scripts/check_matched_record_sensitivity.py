#!/usr/bin/env python3
"""Independently check all sensitivity ranges with dataframe group reductions."""
import json,hashlib
from pathlib import Path
import pandas as pd,numpy as np
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
root=Path('results/structural_comparisons/matched-record-sensitivity-20260927-v1');r=json.loads((root/'receipt.json').read_text());source=Path('results/structural_comparisons/matched-domain-record-summaries-20260927-v1');sr=json.loads((source/'receipt.json').read_text())
assert sha(source/'receipt.json')==r['source_receipt_sha256'] and sha(source/'record_summary.tsv')==sr['artifacts']['record_summary.tsv']
for name,h in r['artifacts'].items():assert sha(root/name)==h
x=pd.read_csv(source/'record_summary.tsv',sep='\t');y=pd.read_csv(root/'order_weighting_sensitivity.tsv',sep='\t');base=['guide','policy','scenario_id','boundary','mask','cohort','screen'];weights=['record','family_equal','taxon_equal'];metrics=['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference'];tol=r['sign_numeric_tolerance'];counts={};values=0
flip=((x[[f'rmsd_difference_{w}_mean' for w in weights]].min(axis=1)<-tol)&(x[[f'rmsd_difference_{w}_mean' for w in weights]].max(axis=1)>tol));flips=x[base].assign(flip=flip).groupby(base).flip.any()
assert int(flips.sum())==r['groups_with_weighting_sign_flip_within_order']
for w in weights:
 actual=y[y.weighting.eq(w)].set_index(base).sort_index();g=x.groupby(base);assert actual.index.is_unique and len(actual)==len(g)==20736
 assert (actual.matched_records==g.matched_records.first()).all() and (actual.weighting_sign_flip_within_any_order==flips.astype(int)).all()
 for metric in metrics:
  low=g[f'{metric}_{w}_mean'].min();high=g[f'{metric}_{w}_mean'].max()
  for suffix,v in [('min',low),('max',high),('range',high-low)]:assert np.allclose(actual[metric+'_order_'+suffix],v,rtol=1e-12,atol=1e-12,equal_nan=True);values+=len(actual)
 low=g[f'rmsd_difference_{w}_mean'].min();high=g[f'rmsd_difference_{w}_mean'].max();status=pd.Series('mixed_sign_or_numeric_zero',index=low.index);status[low>tol]='positive_all_orders';status[high < -tol]='negative_all_orders';status[(low>=-tol)&(high<=tol)]='within_numeric_zero_all_orders';status[g.matched_records.first().eq(0)]='no_matched_records'
 assert (actual.rmsd_order_sign_status==status).all()
 counts.update({w+':'+k:int(v) for k,v in status.value_counts().items()})
assert counts==r['rmsd_sign_counts'] and len(y)==r['weighting_rows']==62208
result=dict(status='passed_full_matched_record_sensitivity_readback',weighting_rows=len(y),numeric_range_values_checked=values,rmsd_sign_counts=counts,groups_with_weighting_sign_flip_within_order=int(flips.sum()),source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),scope='All range endpoints/widths, source cohort counts, weighting-flip flags and sign categories rebuilt by dataframe reductions; no preferred order or biological inference.')
with Path('metadata/matched_record_sensitivity_readback_20260927.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
