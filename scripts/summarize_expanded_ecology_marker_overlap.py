#!/usr/bin/env python3
"""Measure source-specific common markers; no ecological effect inference."""
import csv
import hashlib
import itertools
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    evidence_path = ROOT / 'metadata/species_ecology_evidence_expanded_20260922.tsv'
    evidence_receipt = ROOT / 'metadata/ecology_evidence_expansion_20260922_receipt.json'
    availability_path = ROOT / 'results/structures/current-marker-availability-20260922-v2/marker_availability.tsv'
    availability_receipt = ROOT / 'metadata/current_marker_availability_receipt.json'
    assert sha(evidence_path) == json.loads(evidence_receipt.read_text())['artifacts'][evidence_path.name]
    assert sha(availability_path) == json.loads(availability_receipt.read_text())['artifacts'][availability_path.name]
    evidence = read(evidence_path)
    assert len({r['taxon_id'] for r in evidence}) == len(evidence)
    taxa = {r['taxon_id'] for r in evidence}
    sets = {t: {s: set() for s in ('sequence', 'esmfold', 'alphafold', 'any_model')} for t in taxa}
    seen = set()
    for row in read(availability_path):
        taxon, marker = row['taxon_id'], row['marker']
        if taxon not in taxa:
            continue
        key = taxon, marker
        if key in seen:
            raise ValueError('Multiple proteins per marker/taxon require explicit copy handling')
        seen.add(key)
        status = row['availability']
        assert status in ('both', 'esmfold_only', 'alphafold_only', 'neither_catalog')
        sets[taxon]['sequence'].add(marker)
        if status in ('both', 'esmfold_only'):
            sets[taxon]['esmfold'].add(marker)
        if status in ('both', 'alphafold_only'):
            sets[taxon]['alphafold'].add(marker)
        if status != 'neither_catalog':
            sets[taxon]['any_model'].add(marker)
    pairs = []
    for a, b in itertools.combinations(evidence, 2):
        x, y = sets[a['taxon_id']], sets[b['taxon_id']]
        af, esm = x['alphafold'] & y['alphafold'], x['esmfold'] & y['esmfold']
        any_model = x['any_model'] & y['any_model']
        pairs.append(dict(taxon_a=a['taxon_id'], species_a=a['species_name'], state_a=a['state'],
                          taxon_b=b['taxon_id'], species_b=b['species_name'], state_b=b['state'],
                          different_published_states=a['state'] != b['state'],
                          common_sequences=len(x['sequence'] & y['sequence']),
                          common_esmfold=len(esm), common_alphafold=len(af),
                          common_either_matched_source=len(af | esm),
                          common_any_model=len(any_model),
                          common_only_mixed_sources=len(any_model - (af | esm)),
                          alphafold_marker_ids=';'.join(sorted(af)), esmfold_marker_ids=';'.join(sorted(esm))))
    grouped = defaultdict(list)
    for row in evidence:
        grouped[row['provisional_transition_group']].append(row)
    groups = []
    for group, members in sorted(grouped.items()):
        counts = {}
        for source in ('sequence', 'esmfold', 'alphafold', 'any_model'):
            common = set.intersection(*(sets[r['taxon_id']][source] for r in members))
            counts['common_' + source] = len(common)
        groups.append(dict(provisional_group=group, taxa=len(members),
                           taxon_ids=';'.join(r['taxon_id'] for r in members),
                           published_states=';'.join(sorted({r['state'] for r in members})), **counts))
    artifacts = {}
    for name, rows in [('ecology_expanded_pair_marker_overlap.tsv', pairs),
                       ('ecology_expanded_group_marker_overlap.tsv', groups)]:
        path = ROOT / 'metadata' / name
        with path.open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
        artifacts[name] = sha(path)
    receipt = dict(status='complete_unqualified_source_specific_marker_overlap', taxa=len(taxa),
                   pairwise_comparisons=len(pairs), provisional_groups=len(groups),
                   inputs={str(p.relative_to(ROOT)): sha(p) for p in
                           (evidence_path, evidence_receipt, availability_path, availability_receipt)},
                   script_sha256=sha(Path(__file__)), artifacts=artifacts,
                   scope='Descriptive availability, not confidence-qualified sites, orthology validation, independent transitions, statistical power or ecological effects. Mixed-source-only coverage cannot replace matched-source comparisons.')
    (ROOT / 'metadata/ecology_expanded_marker_overlap_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
