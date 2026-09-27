#!/usr/bin/env python3
"""Independently rejoin selected targets and audit every concentration table row."""
import argparse,hashlib,json,math
from pathlib import Path
import pandas as pd

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
plan=json.loads(a.plan.read_text());root=Path(plan['output']);receipt=root/'receipt.json';r=json.loads(receipt.read_text())
assert r['plan_sha256']==sha(a.plan)
for path,h in plan['pins'].items():assert sha(path)==h
for name,h in r['artifacts'].items():assert sha(root/name)==h
nodes=[]
with open(plan['graph']+'/target_nodes.jsonl') as f:
    for line in f:
        x=json.loads(line);nodes.append({k:x[k] for k in ['node_id','guide','taxon_id','family']})
nodes=pd.DataFrame(nodes)
selected=pd.read_csv(plan['selection']+'/selections.tsv.gz',sep='\t',usecols=['target_id','policy','scenario_id'],dtype=str)
assert not selected.duplicated(['target_id','policy','scenario_id']).any()
joined=selected.merge(nodes,left_on='target_id',right_on='node_id',validate='many_to_one',how='left',indicator=True)
assert (joined['_merge']=='both').all() and len(joined)==r['selected_records']
keys=['guide','policy','scenario_id'];frames=[]
for dimension,col in [('taxon','taxon_id'),('family','family')]:
    actual=joined.groupby(keys+[col],sort=True).size().rename('selected_targets').reset_index()
    base=nodes.groupby(['guide',col]).size().rename('original_targets').reset_index()
    actual=actual.merge(base,on=['guide',col],validate='many_to_one').rename(columns={col:'group_id'})
    actual['dimension']=dimension;frames.append(actual)
expected=pd.concat(frames,ignore_index=True).set_index(keys+['dimension','group_id']).sort_index()
observed=pd.read_csv(root/'group_counts.tsv',sep='\t').set_index(keys+['dimension','group_id']).sort_index()
assert expected.index.equals(observed.index) and len(expected)==r['group_rows']
for col in ['selected_targets','original_targets']:assert (expected[col]==observed[col]).all()
summary=pd.read_csv(root/'concentration.tsv',sep='\t').set_index(keys+['dimension']).sort_index()
coverage=pd.read_csv(plan['balance']+'/selection_coverage.tsv',sep='\t').set_index(keys)
expected_index={(*k,d) for k in coverage.index for d in ['taxon','family']}
assert set(summary.index)==expected_index and len(summary)==2*len(coverage)==r['summary_rows']
groups={k:v.droplevel(keys+['dimension']) for k,v in expected.groupby(level=keys+['dimension'])}
for k,row in summary.iterrows():
    dim=k[3];col='taxon_id' if dim=='taxon' else 'family';subset=nodes[nodes.guide==k[0]]
    group=groups.get(k,expected.iloc[:0])
    counts=group.selected_targets;total=int(counts.sum());original=len(subset);ng=subset[col].nunique()
    values={'original_targets':original,'selected_targets':total,'original_groups':ng,'selected_groups':len(group),'groups_without_matches':ng-len(group),'maximum_group_count':int(counts.max()) if len(counts) else 0,'top_five_share':float(counts.nlargest(5).sum()/total) if total else math.nan,'inverse_concentration':float(total**2/(counts**2).sum()) if total else math.nan}
    for field,value in values.items():assert (math.isnan(value) and pd.isna(row[field])) or math.isclose(row[field],value,rel_tol=1e-12,abs_tol=1e-12),(k,field)
    if total:
        obs=observed.xs(k)
        for field,val in [('share_of_selected',counts/total),('fraction_of_original_group_selected',counts/group.original_targets)]:
            assert all(math.isclose(x,y,rel_tol=1e-12,abs_tol=1e-12) for x,y in zip(obs[field],val)),(k,field)
for path,h in plan['pins'].items():assert sha(path)==h
for name,h in r['artifacts'].items():assert sha(root/name)==h
proof={'status':'passed_full_selected_control_concentration_readback','producer_receipt_sha256':sha(receipt),'selected_records':len(joined),'group_rows':len(expected),'summary_rows':len(summary),'scope':'Independent pandas join/groupby reconstruction of all selected identities, group counts, baseline counts, selection fractions, summary concentration and absent-group counts. No producer helper imported. Descriptive concentration does not establish independent replication or phylogenetic correction.'}
a.output.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof,indent=2))
