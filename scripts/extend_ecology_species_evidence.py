#!/usr/bin/env python3
"""Build a separate expanded evidence table without changing frozen inputs."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    paths = {key: ROOT / value for key, value in {
        'additions': 'config/ecology_species_evidence_additions_20260922.json',
        'previous': 'metadata/species_ecology_evidence.tsv',
        'previous_receipt': 'metadata/species_ecology_evidence_receipt.json',
        'manifest': 'metadata/analysis_manifest.tsv',
        'candidates': 'metadata/ecology_candidates.tsv',
        'coverage': 'metadata/current_marker_taxon_availability.tsv',
    }.items()}
    old_receipt = json.loads(paths['previous_receipt'].read_text())
    assert sha(paths['previous']) == old_receipt['table_sha256']
    assert sha(paths['manifest']) == old_receipt['manifest_sha256']
    assert sha(paths['candidates']) == old_receipt['genus_candidates_sha256']
    evidence = rows(paths['previous'])
    fields = list(evidence[0])
    names = {r['species_name'] for r in evidence}
    config = json.loads(paths['additions'].read_text())
    manifest = rows(paths['manifest'])
    genus = {r['taxon_id']: r for r in rows(paths['candidates'])}
    coverage = {r['taxon_id']: r for r in rows(paths['coverage'])}
    new_coverage = []
    for record in config['records']:
        name = record['species_name']
        assert name not in names, name
        names.add(name)
        matches = [r for r in manifest if r['species_name'] == name and r['study_role'] == 'ingroup']
        assert len(matches) == 1, name
        taxon = matches[0]
        candidate = genus[taxon['taxon_id']]['candidate_primary_lifestyle']
        row = dict.fromkeys(fields, '')
        row.update(record, taxon_id=taxon['taxon_id'], assembly_accession=taxon['assembly_accession'],
                   trait=evidence[0]['trait'], source_species_name=name,
                   evidence='Explicit species classification in the cited Results passage.',
                   evidence_level='published_species_classification',
                   genus_candidate_primary_lifestyle=candidate,
                   conflicts_with_genus_ECM_candidate=str(candidate == 'ectomycorrhizal' and record['state'] == 'asymbiotic'),
                   selected_isolate_experimentally_verified='False',
                   confirmatory_test_status='pending_taxonomy_phylogeny_and_replicated_transition_review')
        for field in ('source_doi', 'source_url', 'source_locator', 'scope_note'):
            row[field] = config[field]
        evidence.append(row)
        c = coverage[taxon['taxon_id']]
        assert c['species_name'] == name
        new_coverage.append(dict(c, curated_state=record['state']))
    artifacts = {}
    for filename, data in [('species_ecology_evidence_expanded_20260922.tsv', evidence),
                           ('ecology_additions_marker_coverage_20260922.tsv', new_coverage)]:
        path = ROOT / 'metadata' / filename
        with path.open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(data)
        artifacts[filename] = sha(path)
    receipt = dict(status='complete_expanded_published_species_classifications',
                   previous_statements=len(evidence)-len(new_coverage), added_statements=len(new_coverage),
                   statements=len(evidence), inputs={k: sha(v) for k, v in paths.items()},
                   script_sha256=sha(Path(__file__)), artifacts=artifacts,
                   scope='Availability precedes confidence filtering. Neither group labels nor tip counts establish independent transitions. Original frozen evidence and downstream analyses are unchanged.')
    (ROOT / 'metadata/ecology_evidence_expansion_20260922_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
