#!/usr/bin/env python3
"""Join full realignment sensitivity to historical flags and repeat the DNA screen.

Historical model flags remain historical: they are neither cleared nor asserted
to have been measured on the new alignment. No selection eligibility is assigned.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

from Bio import SeqIO


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path):
    with Path(path).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    result = {r['case_id']: r for r in rows}
    if len(result) != len(rows):
        raise ValueError(f'Duplicate case: {path}')
    return result


def information(path):
    records = list(SeqIO.parse(path, 'fasta'))
    seqs = [str(r.seq) for r in records]
    ids = [r.id for r in records]
    if not seqs or len(set(ids)) != len(ids) or len({len(s) for s in seqs}) != 1:
        raise ValueError(f'Invalid alignment grid: {path}')
    if set(''.join(seqs)) - set('ACGT?') or len(seqs[0]) % 3:
        raise ValueError(f'Unexpected codon alphabet/dimensions: {path}')
    counts = [Counter(b for b in col if b in 'ACGT') for col in zip(*seqs)]
    informative = sum(sum(n >= 2 for n in c.values()) >= 2 for c in counts)
    distinct = len(set(seqs))
    return ids, dict(taxa=len(seqs), nucleotide_columns=len(seqs[0]),
                     distinct_aligned_sequences=distinct,
                     variable_nucleotide_columns=sum(len(c) > 1 for c in counts),
                     parsimony_informative_nucleotide_columns=informative,
                     information_screen_passes=distinct >= 4 and informative > 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    old = Path('results/cds/genus-codon-diagnostic-inputs-v1')
    new = Path('results/cds/full-group-codon-realignment-projection-20260927-v1')
    sensitivity = Path('results/cds/full-group-alignment-sensitivity-summary-20260927-v1')
    paths = [old / 'receipt.json', new / 'receipt.json', sensitivity / 'receipt.json',
             Path('metadata/full_codon_realignment_projection_readback_20260927.json'),
             Path('metadata/codon_analysis_readiness_receipt.json'),
             Path('metadata/codon_analysis_readiness.tsv'),
             Path('metadata/genus_codon_tree_information.tsv'), sensitivity / 'cases.tsv']
    hashes = {str(p): sha(p) for p in paths}
    old_receipt, new_receipt, sens_receipt, proof, readiness_receipt = [
        json.loads(p.read_text()) for p in paths[:5]]
    if (proof['status'] != 'passed_full_codon_projection_and_correspondence_readback'
            or proof['producer_receipt_sha256'] != sha(new / 'receipt.json')
            or sens_receipt['source_receipt_sha256'] != sha(new / 'receipt.json')
            or sens_receipt['source_readback_sha256'] != sha(paths[3])):
        raise ValueError('Projection proof/source mismatch')
    for p, expected in [(paths[5], readiness_receipt['artifacts']['case_readiness.tsv']),
                        (paths[6], readiness_receipt['sources'][str(paths[6])]),
                        (paths[7], sens_receipt['artifacts']['cases.tsv'])]:
        if sha(p) != expected:
            raise ValueError(f'Source hash mismatch: {p}')
    readiness, original_screen, sensitivities = [table(p) for p in paths[5:]]
    projected = table(new / 'cases.tsv')
    if sha(new / 'cases.tsv') != new_receipt['artifacts']['cases.tsv']:
        raise ValueError('Projection case table changed')
    if not (set(readiness) == set(original_screen) == set(sensitivities) == set(projected)):
        raise ValueError('Case grids differ')
    rows = []
    transitions = Counter()
    for case, previous in readiness.items():
        metrics = []
        names = []
        for root, receipt in [(old, old_receipt), (new, new_receipt)]:
            relative = f'{case}/codons.fna'
            if sha(root / relative) != receipt['artifacts'][relative]:
                raise ValueError(f'Codon input changed: {root / relative}')
            ids, values = information(root / relative)
            metrics.append(values)
            names.append(set(ids))
        if names[0] != names[1]:
            raise ValueError(f'Taxon grid changed: {case}')
        before, after = metrics
        for key, value in before.items():
            if key != 'information_screen_passes' and value != int(original_screen[case][key]):
                raise ValueError(f'Original screen not reproduced: {case}, {key}')
        for key in ['taxa', 'nucleotide_columns', 'distinct_aligned_sequences']:
            if before[key] != int(previous[key]):
                raise ValueError(f'Readiness mismatch: {case}, {key}')
        if before['parsimony_informative_nucleotide_columns'] != int(previous['informative_nucleotide_columns']):
            raise ValueError(f'Readiness information mismatch: {case}')
        if before['information_screen_passes'] != (original_screen[case]['status'] == 'ready_for_supported_tree_diagnostic'):
            raise ValueError(f'Original disposition differs: {case}')
        sens = sensitivities[case]
        if (before['nucleotide_columns'] != 3 * int(sens['original_columns'])
                or after['nucleotide_columns'] != 3 * int(sens['local_retained_columns'])
                or after['taxa'] != int(projected[case]['taxa'])
                or original_screen[case]['marker_copy_caveat'] != previous['copy_caveat']
                or projected[case]['marker_copy_caveat'] != previous['copy_caveat']):
            raise ValueError(f'Sensitivity/copy join mismatch: {case}')
        if previous['selection_eligibility'] != 'not_established':
            raise ValueError(f'Unexpected eligibility: {case}')
        transition = ('pass' if before['information_screen_passes'] else 'fail') + '_to_' + ('pass' if after['information_screen_passes'] else 'fail')
        row = dict(previous)
        row.update({'local_' + k: v for k, v in after.items()})
        row.update({'alignment_' + k: v for k, v in sens.items() if k != 'case_id'})
        row['information_screen_transition'] = transition
        row['historical_flags_scope'] = 'original_alignment_diagnostics_preserved;local_models_not_yet_fitted'
        rows.append(row)
        transitions[transition] += 1
    for p in paths:
        if sha(p) != hashes[str(p)]:
            raise ValueError(f'Input changed during run: {p}')
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / 'cases.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    receipt = dict(status='complete_codon_alignment_information_and_flag_integration',
                   cases=len(rows), information_screen_transitions=dict(transitions),
                   cases_with_preserved_historical_flags=sum(bool(r['case_specific_review_flags']) for r in rows),
                   selection_eligible_cases=0, sources=hashes, script_sha256=sha(__file__),
                   artifacts={'cases.tsv': sha(args.output / 'cases.tsv')},
                   scope='All original flags and fields preserved. Canonical ACGT column counts; distinct aligned strings include missing symbols. Four distinct strings and at least one informative column define the existing information screen, not selection eligibility. No retention cutoff or historical-flag clearance. Local tree/model diagnostics remain required.')
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ['cases', 'information_screen_transitions', 'cases_with_preserved_historical_flags']}))


if __name__ == '__main__':
    main()
