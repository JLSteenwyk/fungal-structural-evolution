#!/usr/bin/env python3
"""Combine unchanged and FCS-omitted normalized exposure observations.

Require unchanged site coordinates; prove affected rows exactly equal their
already-audited sensitivity projection. Does not merge or infer site rates.
"""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table


def signature(row, fields):
    return hashlib.sha256(json.dumps([row[f] for f in fields]).encode()).digest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['baseline', 'changed', 'baseline-inputs', 'inputs', 'changed-readback', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    a = parser.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    receipts = {k: checked_receipt(getattr(a, k)) for k in ['baseline', 'changed', 'baseline_inputs', 'inputs']}
    br, cr, ir = (receipts[k] for k in ['baseline', 'changed', 'inputs'])
    if any(r['status'] != 'complete_reference_normalization' for r in [br, cr]):
        raise ValueError('Incomplete normalization')
    if any(br[k] != cr[k] for k in ['snapshot_receipt_sha256', 'normalization_config_sha256']):
        raise ValueError('Normalization/snapshot differs')
    rb = json.loads(a.changed_readback.read_text())
    if rb['status'] != 'passed_full_accessibility_normalization_readback' or rb['normalized_receipt_sha256'] != sha(a.changed / 'receipt.json'):
        raise ValueError('Matching full changed-row readback required')
    if ir['source_input_receipt_sha256'] != sha(a.baseline_inputs / 'receipt.json'):
        raise ValueError('Input lineage differs')
    disposition = read_table(a.inputs / 'marker_disposition.tsv')
    ready = {r['marker']: r for r in disposition if r['baseline_status'] == 'ready_for_inference'}
    changed = {m for m, r in ready.items() if r['sensitivity_status'] == 'ready_for_inference'}
    unchanged = {m for m, r in ready.items() if r['sensitivity_status'] == 'unchanged_reuse_baseline_fit'}
    if changed | unchanged != set(ready) or any(int(r['all_missing_columns_removed']) for r in ready.values()):
        raise ValueError('Unsupported eligibility or coordinate change')
    omitted = {(r['marker'], r['taxon_id']) for r in read_table(a.inputs / 'omitted_marker_observations.tsv')}
    expected, seen = {}, {}
    columns = {}
    for marker in ready:
        folder = (a.inputs if marker in changed else a.baseline_inputs) / marker
        columns[marker] = [r['matrix_column_1based'] for r in read_table(folder / 'columns.tsv')]
        for record in SeqIO.parse(folder / 'aa.faa', 'fasta'):
            key = marker, record.id
            if key in expected or key in omitted:
                raise ValueError('Duplicate or omitted taxon retained')
            expected[key] = str(record.seq)
            seen[key] = bytearray(len(record.seq))
    partial = {}
    with gzip.open(a.changed / 'normalized_paired_sites.tsv.gz', 'rt') as handle:
        reader = csv.DictReader(handle, delimiter='\t'); fields = reader.fieldnames
        for row in reader:
            key = row['marker'], row['taxon_id'], row['paired_column_1based']
            if key in partial or key[0] not in changed:
                raise ValueError('Invalid changed-row identity')
            partial[key] = signature(row, fields)
    if len(partial) != cr['rows']:
        raise ValueError('Changed row count differs')
    a.output.mkdir(parents=True)
    counts, removed, source_count = Counter(), Counter(), 0
    with gzip.open(a.baseline / 'normalized_paired_sites.tsv.gz', 'rt') as src, gzip.open(a.output / 'normalized_paired_sites.tsv.gz', 'wt') as dst:
        reader = csv.DictReader(src, delimiter='\t')
        if reader.fieldnames != fields:
            raise ValueError('Schemas differ')
        writer = csv.DictWriter(dst, fields, delimiter='\t', lineterminator='\n'); writer.writeheader()
        for row in reader:
            source_count += 1
            key = row['marker'], row['taxon_id']
            if key in omitted:
                removed[key] += 1
                continue
            index = int(row['paired_column_1based']) - 1
            if key not in expected or not 0 <= index < len(expected[key]) or seen[key][index]:
                raise ValueError('Invalid or duplicate retained row')
            if expected[key][index] == '?' or expected[key][index] != row['amino_acid'] or columns[key[0]][index] != row['matrix_column_1based']:
                raise ValueError('Row does not match retained alignment')
            if key[0] in changed and partial.pop((*key, row['paired_column_1based']), None) != signature(row, fields):
                raise ValueError('Changed row differs from audited projection')
            seen[key][index] = 1
            writer.writerow(row); counts[key[0]] += 1
    if partial or source_count != br['rows'] or set(removed) != omitted:
        raise ValueError('Incomplete source, changed or omission grid')
    for key, sequence in expected.items():
        if seen[key] != bytearray(c != '?' for c in sequence):
            raise ValueError('Retained alignment grid incomplete')
    result = {'status': 'complete_full_cohort_fcs_normalized_accessibility',
              'markers': len(ready), 'changed_markers': len(changed), 'unchanged_markers': len(unchanged),
              'rows': sum(counts.values()), 'removed_rows': sum(removed.values()),
              'changed_rows_checked_all_fields': cr['rows'], 'marker_row_counts': dict(sorted(counts.items())),
              'source_receipts': {k: sha(getattr(a, k) / 'receipt.json') for k in receipts},
              'source_paths': {k: str(getattr(a, k)) for k in receipts},
              'changed_readback_sha256': sha(a.changed_readback), 'script_sha256': sha(Path(__file__)),
              'artifacts': {p.name: sha(p) for p in a.output.iterdir()},
              'interpretation': 'Full retained AA observation grid verified; all changed rows equal audited sensitivity normalization and unchanged rows retain baseline values. Whole marker/taxon FCS omission sensitivity, not confirmed contamination removal. No rates or coupling fits merged or completed. Source cohorts overlap and are not independent replicates.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['marker_row_counts','artifacts','source_receipts','source_paths']}, indent=2))


if __name__ == '__main__':
    main()
