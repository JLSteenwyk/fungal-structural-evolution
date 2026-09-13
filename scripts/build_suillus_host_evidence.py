#!/usr/bin/env python3
"""Bind curated host evidence to every sampled Suillus entry without imputation."""
import csv,json
from collections import Counter
from pathlib import Path
from import_ecology_candidates import ROOT,sha


def main():
    config=ROOT/'config/suillus_host_evidence.json';r=json.loads(config.read_text());manifest=ROOT/'metadata/analysis_manifest.tsv'
    with manifest.open() as f:taxa={x['species_name']:x for x in csv.DictReader(f,delimiter='\t') if x['study_role']=='ingroup' and x['species_name'].startswith('Suillus ')}
    names=[x['species_name'] for x in r['records']]
    if len(names)!=len(set(names)) or set(names)!=set(taxa):raise ValueError('Host curation must cover each sampled Suillus entry exactly once')
    out=[]
    for row in r['records']:
        taxon=taxa[row['species_name']];unknown=row['reported_host_groups']=='unknown'
        if unknown!=(row['state_status']!='reported_species_host_classification'):raise ValueError('Unknown host state/status disagree')
        out.append({'taxon_id':taxon['taxon_id'],'assembly_accession':taxon['assembly_accession'],**row,'selected_isolate_experimentally_verified':False,'confirmatory_test_status':'pending_isolate_taxonomy_phylogeny_and_transition_replication'})
    table=ROOT/'metadata/suillus_host_evidence.tsv'
    with table.open('w') as f:
        w=csv.DictWriter(f,list(out[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
    result={'status':'complete_curated_suillus_host_evidence_binding','config_sha256':sha(config),'manifest_sha256':sha(manifest),'script_sha256':sha(Path(__file__)),'taxa':len(out),'state_status_counts':dict(Counter(x['state_status'] for x in out)),'reported_host_group_counts':dict(Counter(x['reported_host_groups'] for x in out)),'table_sha256':sha(table),'interpretation':'Source-reported host classification, not colonization verification of selected isolates. Unknown, multiple hosts and unresolved name-to-table mapping remain explicit. No independent-transition count or host-association test inferred.'}
    (ROOT/'metadata/suillus_host_evidence_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
