#!/usr/bin/env python3
"""Summarize residual gaps from checksum-bound characterization and cache refresh."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gaps', type=Path, required=True)
    parser.add_argument('--refresh', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sources = {}

    def read_table(root, name):
        receipt_bytes = (root / 'receipt.json').read_bytes()
        receipt = json.loads(receipt_bytes)
        raw = (root / name).read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != receipt['artifacts'][name]:
            raise ValueError('Artifact checksum mismatch: ' + name)
        sources[str(root / name)] = digest
        sources[str(root / 'receipt.json')] = hashlib.sha256(receipt_bytes).hexdigest()
        return list(csv.DictReader(raw.decode().splitlines(), delimiter='\t'))

    old = read_table(args.gaps, 'gap_unique_sequences.tsv')
    links = read_table(args.gaps, 'gap_marker_records.tsv')
    refreshed = read_table(args.refresh, 'gap_sequence_availability.tsv')
    old_by = {r['sequence_sha256']: r for r in old}
    new_by = {r['sequence_sha256']: r for r in refreshed}
    assert len(old_by) == len(old) and len(new_by) == len(refreshed)
    assert set(old_by) == set(new_by)
    assert {r['sequence_sha256'] for r in links} == set(old_by)
    for digest, row in new_by.items():
        assert int(row['length']) == int(old_by[digest]['length'])
        assert int(row['verified_cached_candidates']) >= 0
    residual = [r for r in old if int(new_by[r['sequence_sha256']]['verified_cached_candidates']) == 0]
    hashes = {r['sequence_sha256'] for r in residual}
    remaining_links = [r for r in links if r['sequence_sha256'] in hashes]
    categories = Counter()
    for row in residual:
        length = 'above_1024' if int(row['length']) > 1024 else 'at_most_1024'
        alphabet = 'noncanonical' if row['noncanonical_symbols'] else 'canonical'
        categories[length + '_' + alphabet] += 1
    result = dict(
        status='complete_refreshed_marker_gap_summary',
        residual_unique_sequences=len(residual), residual_marker_records=len(remaining_links),
        residual_taxa=len({r['taxon_id'] for r in remaining_links}),
        residual_markers=len({r['marker'] for r in remaining_links}),
        disjoint_length_alphabet_categories=dict(sorted(categories.items())),
        disjoint_length_bins=dict(sorted(Counter(r['length_bin'] for r in residual).items())),
        maximum_length=max(int(r['length']) for r in residual),
        source_sha256=sources, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Residual recovered-marker sequence gaps after the specified frozen cache refresh. Length and alphabet derive from the checksum-bound prior characterization. Not proof of database absence, GPU feasibility, whole-proteome coverage, or authorization to launch predictions.')
    for path, digest in sources.items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
