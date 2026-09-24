#!/usr/bin/env python3
"""Export full candidate domain annotations for focal tree tips and sisters."""
import argparse
import hashlib
import json
from pathlib import Path
import sqlite3


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--context', type=Path, required=True)
    ap.add_argument('--architectures', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    pins = {str(p): sha(p) for p in [a.context, a.architectures, Path(__file__)]}
    context = json.loads(a.context.read_text())
    assert context['status'] == 'complete_focal_resolved_tree_membership_context'
    labels = set()
    for result in context['results']:
        for focal in result['focal_context']:
            labels.add(focal['label'])
            labels.update(focal['sister_labels'])
    con = sqlite3.connect(a.architectures.resolve().as_uri() + '?mode=ro', uri=True)
    con.row_factory = sqlite3.Row
    records = []
    for label in sorted(labels):
        taxon, protein = label.split('_', 1)
        rows = con.execute('SELECT * FROM proteins JOIN queries USING(sequence_id) '
                           'WHERE taxon_id=? AND protein_id=?', (taxon, protein)).fetchall()
        assert len(rows) == 1, (label, len(rows))
        row = dict(rows[0])
        row['candidate_architectures'] = json.loads(row.pop('candidate_architectures_json'))
        row['tree_label'] = label
        records.append(row)
    con.close()
    assert all(sha(p) == digest for p, digest in pins.items())
    result = dict(status='complete_focal_and_sister_candidate_architecture_export',
                  source_hashes=pins, proteins=len(records), records=records,
                  scope='Full retained candidate annotations under all recorded overlap policies. '
                  'Missing retained hits do not establish domain absence, gain, or loss. '
                  'Does not establish function, orthology, or structural acceleration.')
    with a.output.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({'proteins': len(records), 'output': str(a.output)}))


if __name__ == '__main__':
    main()
