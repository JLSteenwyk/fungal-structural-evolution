#!/usr/bin/env python3
"""Select changed paired marker inputs for refitting, preserving original bytes."""
import argparse
import csv
import json
import shutil
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['previous', 'expanded', 'readback', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args(); old = checked_receipt(a.previous); new = checked_receipt(a.expanded)
    audit = json.loads(a.readback.read_text())
    if audit['status'] != 'passed_complete_paired_inputs_from_qualified_arrays_readback' or audit['source_receipt_sha256'] != sha(a.expanded / 'receipt.json'):
        raise ValueError('Expanded input readback incomplete')
    if old['mask'] != new['mask'] or old['eligibility'] != new['eligibility']:
        raise ValueError('Changed filtering')
    def table(root, filename):
        with (root / filename).open() as handle:
            return list(csv.DictReader(handle, delimiter='\t'))
    before = {r['marker']: r for r in table(a.previous, 'marker_summary.tsv')}
    after = {r['marker']: r for r in table(a.expanded, 'marker_summary.tsv')}
    if before.keys() != after.keys():
        raise ValueError('Changed marker universe')
    changed = []
    for marker, row in after.items():
        if row['status'] != 'ready_for_inference':
            continue
        if before[marker]['status'] != 'ready_for_inference' or any(sha(a.previous / marker / name) != sha(a.expanded / marker / name) for name in ['aa.faa', '3di.faa', 'columns.tsv']):
            changed.append(marker)
    if not changed:
        raise ValueError('No changed markers')
    a.output.mkdir(parents=True, exist_ok=False)
    for marker in changed:
        shutil.copytree(a.expanded / marker, a.output / marker)
    for filename in ['marker_summary.tsv', 'taxon_coverage.tsv']:
        rows = [r for r in table(a.expanded, filename) if r['marker'] in changed]
        with (a.output / filename).open('w') as handle:
            w = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
    result = dict(status='complete_changed_marker_input_selection', ready_markers=len(changed), markers=changed,
                  previous_receipt_sha256=sha(a.previous / 'receipt.json'), expanded_receipt_sha256=sha(a.expanded / 'receipt.json'),
                  expanded_readback_sha256=sha(a.readback), script_sha256=sha(Path(__file__)),
                  artifacts={str(f.relative_to(a.output)): sha(f) for f in a.output.rglob('*') if f.is_file()},
                  scope='Exact expanded paired input files for changed markers only. Unchanged marker fits are neither copied nor claimed complete here; full cohort remains the parent snapshot.')
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print('Selected', len(changed), 'changed markers')


if __name__ == '__main__':
    main()
