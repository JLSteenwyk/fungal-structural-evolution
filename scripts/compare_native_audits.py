#!/usr/bin/env python3
"""Demand exact semantic agreement between two complete native feature audits."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new comparison receipt')
    receipts = [checked_receipt(p) for p in [args.reference, args.candidate]]
    if any(r['status'] != 'complete_native_3di_feature_audit' for r in receipts):
        raise ValueError('Complete PAE-qualified feature audits required')
    for key in ['mapping_receipt_sha256', 'pae_receipt_sha256', 'models', 'totals']:
        if receipts[0][key] != receipts[1][key]:
            raise ValueError(f'Different {key}')
    tables = []
    for path in [args.reference, args.candidate]:
        with (path / 'model_summary.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        table = {(r['model_id'], r['version']): r for r in rows}
        if len(table) != len(rows):
            raise ValueError('Repeated model identity')
        tables.append(table)
    if set(tables[0]) != set(tables[1]) or len(tables[0]) != receipts[0]['models']:
        raise ValueError('Different model sets')
    arrays_checked = 0
    for key in tables[0]:
        left, right = (table[key] for table in tables)
        for field in left:
            if field not in ['encoding_path', 'encoding_sha256'] and left[field] != right[field]:
                raise ValueError(f'Different summary field {field}')
        paths = [ROOT / r['encoding_path'] for r in [left, right]]
        for path, row in zip(paths, [left, right]):
            if sha(path) != row['encoding_sha256']:
                raise ValueError('Changed encoding artifact')
        with np.load(paths[0], allow_pickle=False) as a, np.load(paths[1], allow_pickle=False) as b:
            if set(a.files) != set(b.files):
                raise ValueError('Different array fields')
            for name in a.files:
                equal = np.array_equal(a[name], b[name], equal_nan=True) if a[name].dtype.kind == 'f' else np.array_equal(a[name], b[name])
                if not equal:
                    raise ValueError(f'Different array {key}/{name}')
                arrays_checked += 1
    result = {'status': 'exact_semantic_native_audit_agreement', 'models': len(tables[0]),
        'arrays_checked': arrays_checked, 'totals': receipts[0]['totals'],
        'reference_path': str(args.reference), 'candidate_path': str(args.candidate),
        'reference_receipt_sha256': sha(args.reference / 'receipt.json'),
        'candidate_receipt_sha256': sha(args.candidate / 'receipt.json'), 'script_sha256': sha(Path(__file__)),
        'interpretation': 'All model summaries and every encoding array agree exactly, including NaN placement. Regression agreement validates computational equivalence on this complete historical snapshot, not biological accuracy or expanded-dataset completion.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
