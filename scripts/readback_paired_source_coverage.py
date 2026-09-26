#!/usr/bin/env python3
"""Check predictor coverage against actual paired alignment membership."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

from Bio import SeqIO


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    receipt_path = args.comparison / 'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'complete_predictor_source_coverage_comparison'
    pins = {str(receipt_path): sha(receipt_path)}
    for filename, digest in receipt['artifacts'].items():
        path = args.comparison / filename
        assert sha(path) == digest
        pins[str(path)] = digest
    manifest = Path('metadata/analysis_manifest.tsv')
    assert sha(manifest) == receipt['manifest_sha256']
    pins[str(manifest)] = sha(manifest)
    taxa = {r['taxon_id']: r for r in rows(manifest)}
    membership = []
    marker_sets = []
    for source in ['reference', 'local']:
        root = Path(receipt[source + '_path'])
        source_receipt = root / 'receipt.json'
        assert sha(source_receipt) == receipt[source + '_receipt_sha256']
        pins[str(source_receipt)] = sha(source_receipt)
        summary = root / 'marker_summary.tsv'
        pins[str(summary)] = sha(summary)
        markers = rows(summary)
        marker_sets.append({r['marker'] for r in markers})
        assert len(marker_sets[-1]) == len(markers)
        present = set()
        for row in markers:
            if row['status'] != 'ready_for_inference':
                continue
            pair_ids = []
            for filename in ['aa.faa', '3di.faa']:
                path = root / row['marker'] / filename
                pins[str(path)] = sha(path)
                identifiers = [r.id for r in SeqIO.parse(path, 'fasta')]
                assert len(identifiers) == len(set(identifiers))
                assert set(identifiers) <= taxa.keys()
                pair_ids.append(set(identifiers))
            assert pair_ids[0] == pair_ids[1]
            present.update((taxon, row['marker']) for taxon in pair_ids[0])
        membership.append(present)
    assert marker_sets[0] == marker_sets[1]
    assert len(marker_sets[0]) == receipt['markers']
    actual = rows(args.comparison / 'taxon_source_contributions.tsv')
    assert len(actual) == len(taxa) == receipt['taxa']
    assert {r['taxon_id'] for r in actual} == taxa.keys()
    totals = Counter()
    thresholds = Counter()
    missing = []
    for row in actual:
        taxon = row['taxon_id']
        meta = taxa[taxon]
        assert all(row[k] == meta[k] for k in ['species_name', 'study_role'])
        assert row['lineage_group'] == meta['lineage'].split(';')[0]
        a = {m for m in marker_sets[0] if (taxon, m) in membership[0]}
        b = {m for m in marker_sets[0] if (taxon, m) in membership[1]}
        expected = dict(reference_only=len(a-b), local_only=len(b-a), both=len(a&b),
                        neither=len(marker_sets[0]-(a|b)), reference_usable=len(a),
                        local_usable=len(b), union_usable=len(a|b))
        assert all(int(row[k]) == v for k, v in expected.items())
        totals.update(expected)
        for cutoff in [1, 10, 50, 100]:
            thresholds[str(cutoff)] += len(a | b) >= cutoff
        if not a and not b:
            missing.append({k: row[k] for k in ['taxon_id', 'species_name', 'study_role', 'lineage_group']})
    assert receipt['taxon_marker_counts'] == {k: totals[k] for k in receipt['taxon_marker_counts']}
    assert receipt['taxa_with_either'] == thresholds['1']
    assert receipt['taxa_with_reference'] == len({t for t, m in membership[0]})
    assert receipt['taxa_with_local'] == len({t for t, m in membership[1]})
    assert receipt['newly_represented_taxa'] == len({t for t, m in membership[1]} - {t for t, m in membership[0]})
    result = dict(status='passed_complete_alignment_membership_coverage_readback',
                  taxa=len(taxa), markers=len(marker_sets[0]),
                  union_taxon_marker_cells=totals['union_usable'],
                  taxa_with_at_least_n_union_markers=dict(thresholds),
                  taxa_without_usable_markers=missing, source_sha256=pins,
                  script_sha256=sha(Path(__file__)),
                  scope='Checks every taxon-marker membership against both actual AA and 3Di alignment headers. Does not rederive confidence masks, orthology, or calibrated predictor comparability.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}, indent=2))


if __name__ == '__main__':
    main()
