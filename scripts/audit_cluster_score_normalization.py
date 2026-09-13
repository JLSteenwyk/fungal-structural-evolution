#!/usr/bin/env python3
"""Test pinned Foldseek alignment-normalization denominators against saved CIGARs."""
import argparse
import csv
import json
import re
from pathlib import Path
from audit_busco_gene_copies import sha, read_table
from assess_pae_sensitivity import checked_receipt
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['traces', 'review', 'source', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use immutable output')
    reviewed = checked_receipt(a.review)
    source = json.loads((a.source / 'receipt.json').read_text())
    for row in source['files']:
        if sha(Path(row['path'])) != row['sha256']:
            raise ValueError('Source hash mismatch')
    text = (a.source / 'structureconvertalis.cpp').read_text()
    fragment = text.split('case LocalParameters::OUTFMT_ALNTMSCORE:')[1].split('break;')[0]
    if 'std::min(res.qEndPos - res.qStartPos, res.dbEndPos - res.dbStartPos)' not in fragment:
        raise ValueError('Reviewed normalization expression absent')
    scores = {(r['query'], r['target']): r for r in read_table(a.review / 'directed_score_review.tsv')}
    config = json.loads((a.traces / 'config.json').read_text())
    if config['exact_completion_sha256'] != reviewed['exact_completion_sha256']:
        raise ValueError('Trace and score source mismatch')
    fields = config['command'][config['command'].index('--format-output') + 1].split(',')
    rows, seen = [], set()
    for values in csv.reader((a.traces / 'backtraces.tsv').open(), delimiter='\t'):
        if len(values) != len(fields):
            raise ValueError('Trace column mismatch')
        row = dict(zip(fields, values)); key = row['query'], row['target']
        if key in seen or key not in scores:
            raise ValueError('Trace pair mismatch')
        seen.add(key)
        for field in fields[2:-1]:
            if row[field] != scores[key][field]:
                raise ValueError('Trace and score alignment fields differ')
        runs = re.findall(r'(\d+)([MID])', row['cigar'])
        if ''.join(n + op for n, op in runs) != row['cigar']:
            raise ValueError('Unrecognized trace')
        lengths = {op: sum(int(n) for n, code in runs if code == op) for op in 'MID'}
        qspan = int(row['qend']) - int(row['qstart']) + 1
        tspan = int(row['tend']) - int(row['tstart']) + 1
        if lengths['M'] + lengths['I'] != qspan or lengths['M'] + lengths['D'] != tspan or sum(lengths.values()) != int(row['alnlen']):
            raise ValueError('CIGAR does not reconstruct inclusive endpoints')
        denominator = min(qspan, tspan) - 1
        if denominator <= 0:
            raise ValueError('Nonpositive normalization')
        score = float(scores[key]['alntmscore']); bound = lengths['M'] / denominator
        # Export precision may round the exact sum upward by less than 0.0005.
        if score > bound + .0005:
            raise ValueError('Score exceeds denominator-derived bound')
        rows.append(dict(row, matched_positions=lengths['M'], implemented_denominator=denominator,
                         inclusive_min_span=min(qspan, tspan), score=score,
                         denominator_less_than_matched=denominator < lengths['M'],
                         score_above_one=score > 1, theoretical_upper_bound=bound))
    if seen != set(scores):
        raise ValueError('Incomplete trace grid')
    bad = [r for r in rows if r['score_above_one']]
    if not all(r['denominator_less_than_matched'] for r in bad):
        raise ValueError('Overshoot not explained by denominator bound')
    a.output.mkdir(parents=True)
    write_table(a.output / 'normalization_review.tsv', rows)
    write_table(a.output / 'overshoot_cases.tsv', bad)
    receipt = {'status': 'complete_pinned_normalization_source_and_trace_review',
               'source_revision': source['revision'], 'source_receipt_sha256': sha(a.source / 'receipt.json'),
               'trace_config_sha256': sha(a.traces / 'config.json'), 'trace_table_sha256': sha(a.traces / 'backtraces.tsv'),
               'review_receipt_sha256': sha(a.review / 'receipt.json'), 'script_sha256': sha(Path(__file__)),
               'directed_alignments_checked': len(rows), 'scores_above_one': len(bad),
               'denominator_below_matched_positions': sum(r['denominator_less_than_matched'] for r in rows),
               'all_overshoots_within_implemented_normalization_bound': True,
               'interpretation': 'Pinned source uses min(end-start) while every exported CIGAR reconstructs inclusive endpoints. All overshoots have a denominator smaller than their matched-position count and obey its resulting upper bound. This supports an off-by-one normalization explanation; no patched binary or independent TM optimizer has been run. Do not rescale existing scores: normalization also changes distance scale and optimized rotation. Original exclusions remain until corrected rescoring and threshold sensitivity are performed.',
               'artifacts': {f.name: sha(f) for f in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
