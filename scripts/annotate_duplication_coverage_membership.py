#!/usr/bin/env python3
"""Separate excluded candidate-manifest taxa from zero-event reconciled taxa."""
import argparse,csv,json
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--coverage',type=Path,required=True);p.add_argument('--aggregate-readback',type=Path,required=True);p.add_argument('--exclusions',type=Path,required=True);p.add_argument('--profile-tree-audit',type=Path,required=True);p.add_argument('--mafft-tree-audit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
r=json.loads((a.coverage/'receipt.json').read_text());audit=json.loads(a.aggregate_readback.read_text())
if audit['status']!='passed_complete_duplication_sampling_aggregate_readback' or audit['producer_receipt_sha256']!=sha(a.coverage/'receipt.json'):raise ValueError('Coverage aggregate audit differs')
path=a.coverage/'taxon_coverage.tsv'
if sha(path)!=r['artifacts'][path.name]:raise ValueError('Changed coverage table')
exclusions={r['taxon_id']:r for r in json.loads(a.exclusions.read_text())}
pins={str(p):sha(p) for p in [a.coverage/'receipt.json',path,a.aggregate_readback,a.exclusions]};species={}
for guide,ap in [('profile',a.profile_tree_audit),('mafft',a.mafft_tree_audit)]:
    treeaudit=json.loads(ap.read_text());sp=[Path(p) for p in treeaudit['input_hashes'] if Path(p).name=='SpeciesIDs.txt']
    if len(sp)!=1 or sha(sp[0])!=treeaudit['input_hashes'][str(sp[0])]:raise ValueError('Species membership audit binding differs')
    taxa=[line.split(': ',1)[1].rsplit('.',1)[0] for line in sp[0].read_text().splitlines()]
    if len(taxa)!=len(set(taxa)):raise ValueError('Repeated species mapping')
    species[guide]=set(taxa);pins[str(ap)]=sha(ap);pins[str(sp[0])]=sha(sp[0])
a.output.mkdir(parents=True,exist_ok=False);counts={}
with path.open() as f,(a.output/'taxon_coverage_with_membership.tsv').open('w') as target:
    reader=csv.DictReader(f,delimiter='\t');w=csv.DictWriter(target,fieldnames=reader.fieldnames+['in_reconciliation','denominator_status','exclusion_reason'],delimiter='\t',lineterminator='\n');w.writeheader();seen=set()
    for row in reader:
        key=row['guide'],row['taxon']
        if key in seen:raise ValueError('Repeated coverage taxon')
        seen.add(key);included=row['taxon'] in species[row['guide']];n=int(row['terminal_singleton_events'])
        if not included and n:raise ValueError('Excluded taxon has events')
        if not included and row['taxon'] not in exclusions:raise ValueError('Unexplained absent taxon')
        status='events_observed' if n else ('no_qualifying_events' if included else 'not_in_reconciliation')
        label=row['guide']+':'+status;counts[label]=counts.get(label,0)+1
        w.writerow(dict(row,in_reconciliation=int(included),denominator_status=status,exclusion_reason=exclusions[row['taxon']]['reason'] if not included else ''))
for guide,taxa in species.items():
    if not {(guide,t) for t in taxa}<=seen:raise ValueError('Missing reconciled taxon')
for path,h in pins.items():
    if sha(path)!=h:raise ValueError('Changed membership source')
result=dict(status='complete_duplication_coverage_reconciliation_membership',source_sha256=pins,rows=len(seen),counts=counts,artifacts={'taxon_coverage_with_membership.tsv':sha(a.output/'taxon_coverage_with_membership.tsv')},scope='Coverage rows explicitly labeled by actual audited reconciliation species lists; excluded candidate-manifest entries distinguished from observed zero-event species. Original counts retained. No biological absence or missingness correction.')
(a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(counts))
