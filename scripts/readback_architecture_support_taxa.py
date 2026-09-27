#!/usr/bin/env python3
"""Independently regroup all support records using pandas and verify every taxon cell."""
import argparse,json
from pathlib import Path
import pandas as pd
from run_ortholog_pair_guide_comparison import sha
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--summary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
rp=a.summary/'receipt.json';r=json.loads(rp.read_text())
for f,h in r['source_bindings'].items():assert sha(f)==h
out=a.summary/'taxon_support.tsv';assert sha(out)==r['artifacts'][out.name]
keys=['guide','background_set','policy','taxon_id'];factors=['1_25','1_5','2_0']
cols=keys+['family','target_architecture_status','target_same_model']+['within_factor_'+x for x in factors]+['focal_within_factor_'+x for x in factors]
d=pd.read_csv(a.source/'architecture_support.tsv',sep='\t',usecols=cols,keep_default_na=False)
g=d.groupby(keys,sort=True);v=g.size().rename('targets').to_frame();v['families']=g.family.nunique();v['same_model_targets']=g.target_same_model.sum()
for status in ['neither_annotated','one_unannotated','different_ordered_annotations','shared_but_nonconservative','shared_conservative_architecture']:
    v[status]=d.assign(flag=d.target_architecture_status.eq(status)).groupby(keys).flag.sum()
v['supported_families_1_5']=d.loc[d.within_factor_1_5.gt(0)].groupby(keys).family.nunique().reindex(v.index,fill_value=0)
for factor in factors:
    for prefix,source in [('supported_','within_factor_'),('focal_supported_','focal_within_factor_')]:
        v[prefix+factor]=d.assign(flag=d[source+factor].gt(0)).groupby(keys).flag.sum()
actual=pd.read_csv(out,sep='\t').set_index(keys).sort_index()
pd.testing.assert_frame_equal(actual,v[actual.columns],check_dtype=False)
assert len(d)==r['source_rows'] and len(v)==r['taxon_strata']
for row in r['summary']:
    x=v.xs(tuple(row[k] for k in keys[:3]));n=int(x.supported_1_5.sum())
    assert row['targets']==int(x.targets.sum()) and row['target_taxa']==len(x) and row['supported_targets_1_5']==n
    assert row['supported_taxa_1_5']==int(x.supported_1_5.gt(0).sum()) and row['focal_supported_taxa_1_5']==int(x.focal_supported_1_5.gt(0).sum())
    ranked=x.reset_index().sort_values(['supported_1_5','taxon_id'],ascending=[False,True]).head(5)
    assert row['top_five_taxa']==ranked.taxon_id.tolist()
    assert row['top_five_supported_target_fraction']==(int(ranked.supported_1_5.sum())/n if n else None)
for f,h in r['source_bindings'].items():assert sha(f)==h
assert sha(out)==r['artifacts'][out.name]
a.output.write_text(json.dumps(dict(status='passed_full_taxon_architecture_support_readback',producer_receipt_sha256=sha(rp),source_rows=len(d),taxon_strata=len(v),scope='All taxon cells, architecture categories, family counts, distance-range and focal support counts independently regrouped; all concentration summaries verified. This audits availability, not matching or effects.'),indent=2)+'\n')
print(a.output.read_text())
