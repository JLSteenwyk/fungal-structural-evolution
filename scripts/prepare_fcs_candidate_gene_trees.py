#!/usr/bin/env python3
"""Prepare alternative-copy review trees under the existing marker coverage rule."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from Bio import SeqIO

AA = set('ACDEFGHIKLMNPQRSTVWY')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fasta(path):
    records = list(SeqIO.parse(path, 'fasta'))
    result = {r.id: str(r.seq) for r in records}
    if len(result) != len(records):
        raise ValueError('Duplicate identifiers')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--audit', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    source = json.loads((a.source / 'receipt.json').read_text())
    audit = json.loads((a.audit / 'receipt.json').read_text())
    assert audit['status'] == 'complete_independent_stockholm_and_projection_readback'
    assert audit['source_receipt_sha256'] == sha(a.source / 'receipt.json')
    for root, receipt in [(a.source, source), (a.audit, audit)]:
        for name, digest in receipt['artifacts'].items():
            assert sha(root / name) == digest
    a.output.mkdir(parents=True)
    all_rows, marker_rows, artifacts = [], [], {}
    for proof in source['marker_proofs']:
        marker = proof['marker']
        folder = a.source / marker
        for name, digest in proof['artifacts'].items():
            assert sha(folder / name) == digest
        original = fasta(folder / 'baseline_plus_candidates.faa')
        candidates = fasta(folder / 'candidate_projected.faa')
        lengths = {len(s) for s in original.values()}
        assert len(lengths) == 1
        length = lengths.pop()
        minimum = max(50, math.ceil(.3 * length))
        included = {}
        for ident, sequence in original.items():
            observed = sum(c in AA for c in sequence)
            retain = observed >= minimum
            if retain:
                included[ident] = sequence
            all_rows.append({'marker': marker, 'sequence_id': ident,
                             'kind': 'alternative_gene_copy' if ident in candidates else 'baseline',
                             'columns': length, 'observed_residues': observed,
                             'minimum_observed_residues': minimum,
                             'disposition': 'coverage_pass' if retain else 'below_existing_marker_coverage_rule'})
        retained_candidates = set(included) & set(candidates)
        ready = len(included) >= 4 and bool(retained_candidates)
        status = 'ready_for_exploratory_copy_tree' if ready else 'no_coverage_passing_candidate'
        row = {'marker': marker, 'status': status, 'columns': length,
               'baseline_retained': len(included) - len(retained_candidates),
               'candidates_retained': len(retained_candidates),
               'candidates_excluded': len(candidates) - len(retained_candidates),
               'minimum_observed_residues': minimum}
        if ready:
            out = a.output / marker
            out.mkdir()
            path = out / 'input.faa'
            path.write_text(''.join('>' + ident + '\n' + seq + '\n' for ident, seq in sorted(included.items())))
            assert fasta(path) == included
            # Full retained-string comparison, with explicit excluded identifier grid.
            assert all(included[k] == original[k] for k in included)
            row['input_sha256'] = sha(path)
            artifacts[str(path.relative_to(a.output))] = sha(path)
        marker_rows.append(row)
    for name, rows in [('sequence_disposition.tsv', all_rows), ('marker_disposition.tsv', marker_rows)]:
        fields = list(dict.fromkeys(k for row in rows for k in row))
        path = a.output / name
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
        artifacts[name] = sha(path)
    result = {'status': 'complete_candidate_copy_tree_preparation',
              'source_receipt_sha256': sha(a.source / 'receipt.json'),
              'audit_receipt_sha256': sha(a.audit / 'receipt.json'),
              'script_sha256': sha(Path(__file__)),
              'coverage_rule': 'at least max(50, ceil(0.30 * retained marker columns)) canonical residues; same rule for baseline and candidate copies',
              'markers_reviewed': len(marker_rows),
              'markers_ready': sum(r['status'] == 'ready_for_exploratory_copy_tree' for r in marker_rows),
              'candidates_retained': sum(r['candidates_retained'] for r in marker_rows),
              'candidates_excluded': sum(r['candidates_excluded'] for r in marker_rows),
              'artifacts': artifacts,
              'interpretation': 'Exploratory gene-copy trees only. Original selected proteins retained if coverage passes; FCS flags remain relevant. No orthology acceptance or baseline replacement. Below-threshold candidates remain recorded for separate review.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'artifacts'}, indent=2))


if __name__ == '__main__':
    main()
