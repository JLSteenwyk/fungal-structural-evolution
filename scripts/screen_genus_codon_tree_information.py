#!/usr/bin/env python3
"""Reproduce the full within-genus nucleotide-tree information screen."""
import argparse
from collections import Counter
import csv
from pathlib import Path

from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import read_table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checked_receipt(args.inputs)
    fields = ['case_id', 'taxa', 'nucleotide_columns',
              'distinct_aligned_sequences', 'variable_nucleotide_columns',
              'parsimony_informative_nucleotide_columns', 'marker_copy_caveat',
              'status']
    rows = []
    for case in read_table(args.inputs / 'case_summary.tsv'):
        if case['status'] != 'ready_for_tree_and_divergence_diagnostics':
            continue
        records = list(SeqIO.parse(args.inputs / case['case_id'] / 'codons.fna', 'fasta'))
        seqs = [str(record.seq) for record in records]
        expected = int(case['retained_codon_columns']) * 3
        if (len(records) != int(case['retained_taxa']) or
                len({r.id for r in records}) != len(records) or
                any(len(s) != expected for s in seqs)):
            raise ValueError('Alignment grid differs: ' + case['case_id'])
        variable = informative = 0
        for column in zip(*seqs):
            counts = Counter(base for base in column if base in 'ACGT')
            variable += len(counts) > 1
            informative += sum(n >= 2 for n in counts.values()) >= 2
        distinct = len(set(seqs))
        ready = distinct >= 4 and informative > 0
        rows.append(dict(zip(fields, [case['case_id'], len(seqs), expected,
            distinct, variable, informative, case['marker_copy_caveat'],
            'ready_for_supported_tree_diagnostic' if ready else
            'insufficient_sequence_information_for_this_tree_workflow'])))
    with args.output.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print(dict(Counter(row['status'] for row in rows)))


if __name__ == '__main__':
    main()
