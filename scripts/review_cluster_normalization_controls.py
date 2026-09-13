#!/usr/bin/env python3
"""Separate toolchain and normalization effects across complete fixed-alignment controls."""
import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha
from prepare_paired_phylogenetic_inputs import write_table

BOUNDED = ['qcov', 'tcov', 'alntmscore', 'qtmscore', 'ttmscore', 'lddt']
SCORES = ['alntmscore', 'qtmscore', 'ttmscore', 'lddt', 'rmsd']


def decision(row):
    nums = {k: float(v) for k, v in row.items() if k not in ['query', 'target']}
    if any(not math.isfinite(v) or v < 0 for v in nums.values()):
        raise ValueError('Nonfinite or negative score')
    valid = all(nums[k] <= 1 for k in BOUNDED)
    return valid, valid and nums['evalue'] <= 1e-5 and min(nums['qcov'], nums['tcov']) >= .8 and nums['alntmscore'] >= .5


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['controls', 'original', 'output']:
        p.add_argument('--' + name, required=True, type=Path)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use immutable output')
    receipt = checked_receipt(a.controls)
    config = json.loads((a.controls / 'config.json').read_text())
    if receipt['config_sha256'] != sha(a.controls / 'config.json'):
        raise ValueError('Control configuration mismatch')
    original_config = json.loads((a.original / 'config.json').read_text())
    completion = json.loads((a.original / 'completion.json').read_text())
    if completion['alignment_table_sha256'] != sha(a.original / 'alignments.tsv') or completion['config_sha256'] != sha(a.original / 'config.json') or config['exact_completion_sha256'] != sha(a.original / 'completion.json'):
        raise ValueError('Original source mismatch')
    command = original_config['command']
    fields = command[command.index('--format-output') + 1].split(',')
    tables = {}
    for name, path in [('production', a.original / 'alignments.tsv'), ('unmodified', a.controls / 'unmodified.tsv'), ('inclusive-span', a.controls / 'inclusive-span.tsv')]:
        rows = {}
        for values in csv.reader(path.open(), delimiter='\t'):
            if len(values) != len(fields):
                raise ValueError('Column mismatch')
            row = dict(zip(fields, values)); key = row['query'], row['target']
            if key in rows:
                raise ValueError('Duplicate pair')
            decision(row)
            rows[key] = row
        tables[name] = rows
    keys = set(tables['production'])
    if not keys or any(set(table) != keys for table in tables.values()):
        raise ValueError('Pair universe mismatch')
    fixed = [f for f in fields if f not in SCORES]
    output, stage_summaries = [], []
    for stage, before_name, after_name in [('toolchain', 'production', 'unmodified'), ('normalization', 'unmodified', 'inclusive-span')]:
        changed = Counter(); maximum = {f: 0.0 for f in SCORES}
        transitions = Counter(); invalid_before = invalid_after = 0
        for key in sorted(keys):
            before, after = tables[before_name][key], tables[after_name][key]
            if any(before[f] != after[f] for f in fixed):
                raise ValueError('Fixed alignment fields changed')
            if stage == 'normalization' and any(before[f] != after[f] for f in SCORES if f != 'alntmscore'):
                raise ValueError('Normalization patch changed another score output')
            bv, bp = decision(before); av, ap = decision(after)
            invalid_before += not bv; invalid_after += not av
            transitions[str(bp) + '_to_' + str(ap)] += 1
            for field in SCORES:
                delta = abs(float(after[field]) - float(before[field]))
                changed[field] += delta != 0
                maximum[field] = max(maximum[field], delta)
            output.append({'stage': stage, 'query': key[0], 'target': key[1],
                           'before_alntmscore': before['alntmscore'], 'after_alntmscore': after['alntmscore'],
                           'before_valid': bv, 'after_valid': av, 'before_passes': bp, 'after_passes': ap})
        stage_summaries.append({'stage': stage, 'directed_pairs': len(keys),
                                'invalid_rows_before': invalid_before, 'invalid_rows_after': invalid_after,
                                'decision_transitions': dict(transitions), 'score_rows_changed': dict(changed),
                                'maximum_absolute_score_changes': maximum})
    a.output.mkdir(parents=True)
    write_table(a.output / 'directed_control_comparisons.tsv', output)
    result = {'status': 'complete_toolchain_and_normalization_control_readback',
              'control_receipt_sha256': sha(a.controls / 'receipt.json'),
              'original_completion_sha256': sha(a.original / 'completion.json'),
              'script_sha256': sha(Path(__file__)), 'stage_summaries': stage_summaries,
              'interpretation': 'Complete fixed-alignment pair grid and finite/nonnegative scores checked. Toolchain effects are separate from the single output-denominator patch; other score fields must be identical between matched builds. Normalization sensitivity is not an independent structural alignment, orthology, prediction accuracy or novelty validation. Derived memberships require a separate explicit review.',
              'artifacts': {f.name: sha(f) for f in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
