#!/usr/bin/env python3
"""Independently check the full joined ledger with pandas and NumPy DNA counts."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frame(path):
    table = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    assert not table.case_id.duplicated().any(), path
    return table.set_index('case_id').sort_index()


def manual_metrics(path):
    sequences = []
    ids = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            ids.append(line[1:].split()[0])
            sequences.append('')
        elif line.strip():
            sequences[-1] += line.strip()
    assert len(set(ids)) == len(ids)
    assert len({len(s) for s in sequences}) == 1
    array = np.array([list(s) for s in sequences])
    counts = np.stack([(array == b).sum(axis=0) for b in 'ACGT'])
    informative = int(((counts >= 2).sum(axis=0) >= 2).sum())
    return dict(taxa=len(ids), nucleotide_columns=array.shape[1],
                distinct_aligned_sequences=len(set(sequences)),
                variable_nucleotide_columns=int(((counts > 0).sum(axis=0) > 1).sum()),
                parsimony_informative_nucleotide_columns=informative,
                information_screen_passes=len(set(sequences)) >= 4 and informative > 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads((args.input / 'receipt.json').read_text())
    assert digest(args.input / 'cases.tsv') == receipt['artifacts']['cases.tsv']
    for name, expected in receipt['sources'].items():
        assert digest(name) == expected, name
    actual = frame(args.input / 'cases.tsv')
    old = frame('metadata/codon_analysis_readiness.tsv')
    sensitivity = frame('results/cds/full-group-alignment-sensitivity-summary-20260927-v1/cases.tsv')
    pd.testing.assert_frame_equal(actual[old.columns], old)
    for col in sensitivity:
        assert actual['alignment_' + col].equals(sensitivity[col]), col
    transitions = Counter()
    for case, row in actual.iterrows():
        values = manual_metrics(Path('results/cds/full-group-codon-realignment-projection-20260927-v1') / case / 'codons.fna')
        for key, value in values.items():
            assert row['local_' + key] == str(value), (case, key)
        previous = int(row.distinct_aligned_sequences) >= 4 and int(row.informative_nucleotide_columns) > 0
        transition = ('pass' if previous else 'fail') + '_to_' + ('pass' if values['information_screen_passes'] else 'fail')
        assert row.information_screen_transition == transition, case
        assert row.selection_eligibility == 'not_established'
        assert row.historical_flags_scope == 'original_alignment_diagnostics_preserved;local_models_not_yet_fitted'
        transitions[transition] += 1
    assert dict(transitions) == receipt['information_screen_transitions']
    flagged = int((actual.case_specific_review_flags != '').sum())
    assert flagged == receipt['cases_with_preserved_historical_flags']
    assert len(actual) == receipt['cases'] == 1712
    proof = dict(status='passed_full_codon_alignment_integration_readback', cases=len(actual),
                 source_receipt_sha256=digest(args.input / 'receipt.json'),
                 script_sha256=digest(__file__), information_screen_transitions=dict(transitions),
                 preserved_flagged_cases=flagged,
                 scope='Every original ledger field and sensitivity field checked by exact pandas joins. Manual FASTA parser and NumPy base counts independently reproduce all local information metrics and transitions. Does not establish homology correctness or selection eligibility.')
    with args.output.open('x') as handle:
        json.dump(proof, handle, indent=2)
        handle.write('\n')
    print(json.dumps(proof))


if __name__ == '__main__':
    main()
