#!/usr/bin/env python3
"""Attach manually reviewed primary-source ecological statements to selected taxa."""
import csv
import json
from pathlib import Path
from import_ecology_candidates import ROOT, sha


def main():
    config = ROOT / 'config/ecology_species_evidence.json'
    evidence = json.loads(config.read_text())
    manifest = ROOT / 'metadata/analysis_manifest.tsv'
    with manifest.open() as handle:
        taxa = list(csv.DictReader(handle, delimiter='\t'))
    candidates = ROOT / 'metadata/ecology_candidates.tsv'
    with candidates.open() as handle:
        genus = {r['taxon_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    receipt = json.loads((ROOT / 'metadata/ecology_candidates_receipt.json').read_text())
    if sha(candidates) != receipt['table_sha256']:
        raise ValueError('Changed genus candidate evidence')
    rows = []
    names = set()
    for record in evidence['records']:
        name = record['species_name']
        if name in names:
            raise ValueError('Duplicate species statement; use explicit multi-source evidence handling')
        names.add(name)
        matches = [r for r in taxa if r['species_name'] == name and r['study_role'] == 'ingroup']
        if len(matches) != 1:
            raise ValueError('Species statement does not match exactly one fungal entry')
        taxon = matches[0]
        candidate = genus[taxon['taxon_id']]['candidate_primary_lifestyle']
        rows.append({'taxon_id': taxon['taxon_id'], 'assembly_accession': taxon['assembly_accession'],
            'trait': evidence['trait'], **record, 'evidence_level': 'published_species_classification',
            'genus_candidate_primary_lifestyle': candidate,
            'conflicts_with_genus_ECM_candidate': candidate == 'ectomycorrhizal' and record['state'] == 'asymbiotic',
            'selected_isolate_experimentally_verified': False,
            'confirmatory_test_status': 'pending_taxonomy_phylogeny_and_replicated_transition_review'})
    path = ROOT / 'metadata/species_ecology_evidence.tsv'
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    result = {'config_sha256': sha(config), 'manifest_sha256': sha(manifest),
        'genus_candidates_sha256': sha(candidates), 'script_sha256': sha(Path(__file__)),
        'statements': len(rows), 'species': len(names),
        'conflicting_genus_ECM_candidates': sum(r['conflicts_with_genus_ECM_candidate'] for r in rows),
        'table_sha256': sha(path),
        'interpretation': 'Manually reviewed published species classifications. This is not completed ecological coverage, exact-isolate experimental validation, or an independent-transition count.'}
    (ROOT / 'metadata/species_ecology_evidence_receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
