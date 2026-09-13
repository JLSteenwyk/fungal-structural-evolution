#!/usr/bin/env python3
"""Flag explicit hybrid and uncertain species labels without inferring taxonomy."""
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    source = ROOT / 'metadata/analysis_manifest.tsv'
    with source.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    names = Counter(r['species_name'] for r in rows)
    output = []
    for row in rows:
        flags = []
        if re.search(r'\s[x×]\s', row['species_name']):
            flags.append('explicit_hybrid_label')
        if re.search(r'\b(?:sp|cf|aff)\.', row['species_name']):
            flags.append('incompletely_identified_species_label')
        if names[row['species_name']] > 1:
            flags.append('repeated_species_label')
        if flags:
            output.append({k: row[k] for k in ['taxon_id', 'species_name', 'study_role', 'assembly_accession']}
                          | {'flags': ';'.join(flags), 'review_status': 'pending',
                             'analysis_policy': 'Retain raw data; review species identity and test taxon exclusion before final species-level inference'})
    path = ROOT / 'metadata/taxon_label_review.tsv'
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, ['taxon_id', 'species_name', 'study_role', 'assembly_accession',
                                       'flags', 'review_status', 'analysis_policy'], delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(output)
    receipt = {'manifest_sha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'taxa_screened': len(rows),
        'flagged_taxa': len(output), 'flag_counts': dict(Counter(f for r in output for f in r['flags'].split(';'))),
        'interpretation': 'Name-based triage only. Unflagged labels do not establish accepted taxonomy, independent species or absence of hybrid ancestry.',
        'artifact_sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    (ROOT / 'metadata/taxon_label_review_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
