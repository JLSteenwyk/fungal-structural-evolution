#!/usr/bin/env python3
"""Check that expanded paired inputs preserve every prior observed character."""
import argparse
import csv
import json
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['previous', 'expanded', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    old_receipt, new_receipt = checked_receipt(a.previous), checked_receipt(a.expanded)
    for field in ['mask', 'eligibility']:
        if old_receipt[field] != new_receipt[field]:
            raise ValueError('Changed filtering')
    if old_receipt['source_receipts']['matrix'] != new_receipt['source_receipts']['matrix']:
        raise ValueError('Different source matrix')
    coverage = [{(r['marker'], r['taxon_id']): r for r in table(root / 'taxon_coverage.tsv')} for root in [a.previous, a.expanded]]
    old, new = coverage
    if old.keys() != new.keys():
        raise ValueError('Different coverage universe')
    changes = []
    for key in sorted(old):
        if old[key] == new[key]:
            continue
        if int(new[key]['observed']) < int(old[key]['observed']) or (old[key]['taxon_eligible'] == 'True' and new[key]['taxon_eligible'] != 'True'):
            raise ValueError('Lost coverage')
        changes.append(dict(marker=key[0], taxon_id=key[1], previous_observed=int(old[key]['observed']), expanded_observed=int(new[key]['observed']), previous_eligible=old[key]['taxon_eligible'], expanded_eligible=new[key]['taxon_eligible']))
    unchanged_markers = changed_markers = checked_characters = 0
    for summary in table(a.previous / 'marker_summary.tsv'):
        if summary['status'] != 'ready_for_inference':
            continue
        marker = summary['marker']; prior = a.previous / marker; expanded = a.expanded / marker
        columns = [[int(r['matrix_column_1based']) for r in table(root / 'columns.tsv')] for root in [prior, expanded]]
        index = {column: i for i, column in enumerate(columns[1])}
        if not set(columns[0]) <= set(index):
            raise ValueError('Lost alignment columns')
        identical = columns[0] == columns[1]
        for filename in ['aa.faa', '3di.faa']:
            seqs = [{r.id: str(r.seq) for r in SeqIO.parse(root / filename, 'fasta')} for root in [prior, expanded]]
            if not seqs[0].keys() <= seqs[1].keys():
                raise ValueError('Lost eligible taxa')
            identical = identical and seqs[0] == seqs[1]
            for taxon, sequence in seqs[0].items():
                for i, char in enumerate(sequence):
                    value = seqs[1][taxon][index[columns[0][i]]]
                    if char != '?' and value != char:
                        raise ValueError('Changed previously observed character')
                    checked_characters += char != '?'
        unchanged_markers += identical; changed_markers += not identical
    result = dict(status='passed_expanded_inputs_preserve_previous_observations',
                  previous_receipt_sha256=sha(a.previous / 'receipt.json'), expanded_receipt_sha256=sha(a.expanded / 'receipt.json'),
                  script_sha256=sha(Path(__file__)), unchanged_markers=unchanged_markers, changed_markers=changed_markers,
                  checked_previous_observed_characters=checked_characters, coverage_changes=changes,
                  additional_observed_paired_cells=sum(r['expanded_observed'] - r['previous_observed'] for r in changes),
                  scope='Same filtering and matrix; all prior observed AA/3Di characters preserved at matching source columns. Coverage changes are not evolutionary effects; full array readback remains a separate check.')
    with a.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'coverage_changes'}))


if __name__ == '__main__':
    main()
