#!/usr/bin/env python3
"""Quantify taxonomic concentration of available within-architecture controls."""
import argparse,csv,json
from collections import Counter,defaultdict
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--proof',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
rp=a.source/'receipt.json';receipt=json.loads(rp.read_text());proof=json.loads(a.proof.read_text());table=a.source/'architecture_support.tsv'
assert proof['status']=='passed_full_architecture_matched_support_readback' and proof['producer_receipt_sha256']==sha(rp)
assert sha(table)==receipt['artifacts'][table.name]
bindings={str(f):sha(f) for f in [rp,a.proof,table]}
groups=defaultdict(Counter);families=defaultdict(set);supported=defaultdict(set);rows=0
fields=['guide','background_set','policy','taxon_id']
statuses=['neither_annotated','one_unannotated','different_ordered_annotations','shared_but_nonconservative','shared_conservative_architecture']
for row in csv.DictReader(table.open(),delimiter='\t'):
    key=tuple(row[k] for k in fields);c=groups[key];rows+=1;c['targets']+=1
    assert row['target_architecture_status'] in statuses
    c[row['target_architecture_status']]+=1;c['same_model_targets']+=int(row['target_same_model'])
    families[key].add(row['family'])
    for factor in ['1_25','1_5','2_0']:
        n=int(row['within_factor_'+factor]);f=int(row['focal_within_factor_'+factor]);o=int(row['nonfocal_within_factor_'+factor]);assert n==f+o
        c['supported_'+factor]+=n>0;c['focal_supported_'+factor]+=f>0
    if int(row['within_factor_1_5']):supported[key].add(row['family'])
assert rows==receipt['rows']==proof['rows']
a.output.mkdir(parents=True,exist_ok=False);out=a.output/'taxon_support.tsv'
metrics=['targets','families','same_model_targets',*statuses,'supported_families_1_5',*[s+v for v in ['1_25','1_5','2_0'] for s in ['supported_','focal_supported_']]]
with out.open('w') as f:
    w=csv.DictWriter(f,fieldnames=fields+metrics,delimiter='\t',lineterminator='\n');w.writeheader()
    for key,c in sorted(groups.items()):
        c['families']=len(families[key]);c['supported_families_1_5']=len(supported[key]);w.writerow(dict(zip(fields,key))|{m:c[m] for m in metrics})
summary=[]
for group in sorted({k[:3] for k in groups}):
    cells=[(k[3],v) for k,v in groups.items() if k[:3]==group];total=sum(v['targets'] for _,v in cells);support=sum(v['supported_1_5'] for _,v in cells)
    ranked=sorted(cells,key=lambda x:(-x[1]['supported_1_5'],x[0]));summary.append(dict(zip(fields[:3],group))|dict(target_taxa=len(cells),targets=total,supported_targets_1_5=support,supported_taxa_1_5=sum(v['supported_1_5']>0 for _,v in cells),focal_supported_taxa_1_5=sum(v['focal_supported_1_5']>0 for _,v in cells),top_five_supported_target_fraction=sum(v['supported_1_5'] for _,v in ranked[:5])/support if support else None,top_five_taxa=[t for t,_ in ranked[:5]]))
for f,h in bindings.items():assert sha(f)==h
result=dict(status='complete_taxon_architecture_support_summary_pending_readback',source_bindings=bindings,source_rows=rows,taxon_strata=len(groups),summary=summary,artifacts={out.name:sha(out)},scope='Availability only: each guide/set/policy has its own denominator of modeled terminal duplicate targets, not all fungal genes or taxa. Zero-support target taxa retained; entirely absent taxa are outside this observed denominator. Guides are sensitivity alternatives, not replicates. Identical-model targets retained explicitly. No actual matching, structural effects, causal inference or clade coverage claim.')
(a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(source_rows=rows,taxon_strata=len(groups))))
