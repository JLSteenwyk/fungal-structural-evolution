#!/usr/bin/env python3
"""Hold matched model-pair membership fixed across all confidence alternatives."""
import argparse,csv,json,statistics
from collections import defaultdict
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha
from prepare_paired_phylogenetic_inputs import write_table

p=argparse.ArgumentParser(description=__doc__)
for name in ['comparison','readback','output']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();rp=a.comparison/'receipt.json';r=json.loads(rp.read_text());audit=json.loads(a.readback.read_text())
if audit['status']!='passed_full_overlap_model_feature_comparison_readback' or audit['producer_receipt_sha256']!=sha(rp):raise ValueError('Unverified comparison')
table=a.comparison/'model_pair_comparisons.tsv'
if sha(table)!=r['artifacts'][table.name]:raise ValueError('Changed comparisons')
groups=defaultdict(dict);settings={(str(p),str(q)) for p in [0,70,90] for q in ['unfiltered',5,10,15]}
for row in csv.DictReader(table.open(),delimiter='\t'):
 key=tuple(row[k] for k in ['reference_model_id','reference_version','local_model_id','local_version']);setting=row['plddt_cutoff'],row['pae_cutoff']
 if setting in groups[key]:raise ValueError('Repeated model/threshold')
 groups[key][setting]=row
if len(groups)!=r['model_pairs'] or any(set(v)!=settings for v in groups.values()):raise ValueError('Incomplete model/threshold scope')
retained={k:v for k,v in groups.items() if all(x['status']=='compared' for x in v.values())}
if not retained:raise ValueError('No common eligible cohort')
a.output.mkdir(parents=True,exist_ok=False)
write_table(a.output/'model_pair_membership.tsv',[dict(reference_model_id=k[0],reference_version=k[1],local_model_id=k[2],local_version=k[3],eligible_all_thresholds=int(k in retained),eligible_thresholds=sum(x['status']=='compared' for x in v.values())) for k,v in sorted(groups.items())])
summary=[]
for cutoff in [0,70,90]:
 for pae in ['unfiltered',5,10,15]:
  selected=[v[str(cutoff),str(pae)] for v in retained.values()];n=sum(int(x['retained_residues']) for x in selected);diff=sum(int(x['state_mismatches']) for x in selected)
  summary.append(dict(plddt_cutoff=cutoff,pae_cutoff=pae,model_pairs=len(selected),retained_residues=n,state_mismatches=diff,pooled_state_mismatch_fraction=diff/n,median_pair_mismatch_fraction=statistics.median(float(x['state_mismatch_fraction']) for x in selected),mean_pair_mismatch_fraction=statistics.mean(float(x['state_mismatch_fraction']) for x in selected)))
write_table(a.output/'fixed_cohort_threshold_summary.tsv',summary)
changes=[]
for k,v in sorted(retained.items()):
 base=v['0','unfiltered'];medium=v['70','10'];strict=v['90','10']
 changes.append(dict(reference_model_id=k[0],reference_version=k[1],local_model_id=k[2],local_version=k[3],valid_only_fraction=base['state_mismatch_fraction'],plddt70_pae10_fraction=medium['state_mismatch_fraction'],plddt90_pae10_fraction=strict['state_mismatch_fraction'],strict_minus_valid_only=float(strict['state_mismatch_fraction'])-float(base['state_mismatch_fraction']),strict_minus_plddt70=float(strict['state_mismatch_fraction'])-float(medium['state_mismatch_fraction'])))
write_table(a.output/'paired_threshold_changes.tsv',changes)
result=dict(status='complete_fixed_cohort_overlap_sensitivity',source_receipt_sha256=sha(rp),readback_sha256=sha(a.readback),source_table_sha256=sha(table),script_sha256=sha(__file__),original_model_pairs=len(groups),fixed_model_pairs=len(retained),excluded_from_fixed_cohort=len(groups)-len(retained),artifacts={x.name:sha(x) for x in a.output.iterdir()},scope='Same model pairs across twelve confidence alternatives; residue sets still change. Highly selected intersection cohort, not a representative sample or causal confidence intervention. Pair-specific changes descriptive; no independence, error calibration or significance claim.')
(a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print('Fixed pairs:',len(retained),'of',len(groups))
for row in summary:
 if row['pae_cutoff']=='unfiltered' and row['plddt_cutoff']==0 or row['pae_cutoff']==10:print(row)
