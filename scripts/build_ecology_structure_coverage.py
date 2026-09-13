#!/usr/bin/env python3
"""Measure predictor-specific structural coverage of curated ecological comparisons."""
import csv,json
from collections import defaultdict
from pathlib import Path
from import_ecology_candidates import ROOT,sha


def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))

def write(p,rows):
    with p.open('w') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def main():
    ecological=ROOT/'metadata/species_ecology_evidence.tsv';hosts=ROOT/'metadata/suillus_host_evidence.tsv'
    for table,receipt in [(ecological,'species_ecology_evidence_receipt.json'),(hosts,'suillus_host_evidence_receipt.json')]:
        if sha(table)!=json.loads((ROOT/'metadata'/receipt).read_text())['table_sha256']:raise ValueError('Changed ecological evidence')
    evidence=read(ecological);host={r['taxon_id']:r for r in read(hosts)};groups=defaultdict(list)
    for r in evidence:groups[r['provisional_transition_group']].append(r)
    totals={};source_receipts={};contrasts=[]
    for source,folder in [('AlphaFold','paired-inputs-gdm-expanded-v1'),('ESMFold','paired-inputs-esmfold-partial-v1')]:
        b=ROOT/'results/phylogeny'/folder;r=json.loads((b/'receipt.json').read_text())
        for n in ['taxon_coverage.tsv','marker_summary.tsv']:
            if sha(b/n)!=r['artifacts'][n]:raise ValueError('Changed paired input coverage')
        source_receipts[source]=sha(b/'receipt.json')
        ready={r['marker'] for r in read(b/'marker_summary.tsv') if r['status']=='ready_for_inference'};eligible=defaultdict(set)
        for r in read(b/'taxon_coverage.tsv'):
            if r['taxon_eligible']=='True' and r['marker'] in ready:eligible[r['taxon_id']].add(r['marker'])
        totals[source]=eligible
        for group,members in groups.items():
            states={r['state'] for r in members}
            if len(states)<2:continue
            common=set.intersection(*(eligible[r['taxon_id']] for r in members))
            contrasts.append({'provisional_group':group,'source':source,'taxa':len(members),'source_states':';'.join(sorted(states)),'markers_eligible_in_every_group_taxon':len(common),'markers':';'.join(sorted(common)),'interpretation':'Intersection of source-specific eligible coverage only; no verified transition, effect or statistical-power claim.'})
    rows=[{'taxon_id':r['taxon_id'],'species_name':r['species_name'],'source_classification':r['state'],'provisional_transition_group':r['provisional_transition_group'],
           'AlphaFold_usable_markers':len(totals['AlphaFold'][r['taxon_id']]),'ESMFold_usable_markers':len(totals['ESMFold'][r['taxon_id']]),
           'reported_host_groups':host.get(r['taxon_id'],{}).get('reported_host_groups','not_curated'),'confirmatory_test_status':r['confirmatory_test_status']} for r in evidence]
    table=ROOT/'metadata/ecology_paired_coverage.tsv';contrast=ROOT/'metadata/ecology_contrast_common_marker_coverage.tsv';write(table,rows);write(contrast,contrasts)
    result={'status':'complete_ecology_paired_coverage_join','ecology_evidence_sha256':sha(ecological),'host_evidence_sha256':sha(hosts),'paired_input_receipt_sha256':source_receipts,'script_sha256':sha(Path(__file__)),
            'ecologically_curated_species':len(rows),'with_any_source_usable_marker':sum(r['AlphaFold_usable_markers']+r['ESMFold_usable_markers']>0 for r in rows),'table_sha256':sha(table),'contrast_table_sha256':sha(contrast),
            'interpretation':'Availability only. Separate-source coverage exposes potential prediction-source confounding. No pooled inference, independent-origin count, ecological effect or adequacy/power claim.'}
    (ROOT/'metadata/ecology_paired_coverage_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
