#!/usr/bin/env python3
"""Known choices, complete synthetic handoff and rehashed export corruption checks."""
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from compare_sequence_first_references import FIELDS, comparison, validate_row
from readback_sequence_first_reference_comparison import expected
from screen_duplication_alignment_reuse import sha


def row(n, ties, model='', status='no_modeled_nonfocal_sister'):
    return dict(family='OG0', gene_a=f'T{n}_a', gene_b=f'T{n}_b', status=status,
                sequence_first_status='fixture_context', sequence_first_nearest_genes=json.dumps(ties),
                sequence_first_chosen_gene=ties[0] if ties else '', sequence_first_chosen_model=model,
                sequence_first_chosen_version='10' if model else '', sequence_first_minimum_distance='1.5',
                availability_distance_increment='0' if model else '')


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def table(path, rows):
    with path.open('w') as handle:
        w = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    cases = [
        (row(1, ['A']), row(1, ['A']), 'same_nearest_sequence_gene_set', '1', '0', ''),
        (row(2, ['A'], 'M', 'provisional_reference_available'),
         row(2, ['A'], 'M', 'provisional_reference_available'), 'same_nearest_sequence_gene_set', '1', '1', '1'),
        (row(3, ['A'], 'M1'), row(3, ['B'], 'M2'), 'disjoint_nearest_sequence_gene_sets', '0', '1', '0'),
        (row(4, ['A', 'B']), row(4, ['B', 'C'], 'M'), 'overlapping_nearest_sequence_gene_sets', '0', '0', ''),
        (row(5, ['A'], 'M'), row(5, ['B'], 'M'), 'disjoint_nearest_sequence_gene_sets', '0', '1', '1'),
        (row(6, ['A'], status='parent_reported_duplication'), row(6, ['A']),
         'ineligible_sequence_context_in_one_or_both', '', '', ''),
        (row(7, ['A']), None, 'profile_only_candidate', '', '', ''),
        (None, row(8, ['A']), 'mafft_only_candidate', '', '', ''),
        (row(9, [], status='no_parent'), row(9, [], status='no_parent'),
         'ineligible_sequence_context_in_one_or_both', '', '', ''),
    ]
    left, right, output = [], [], []
    for a, b, relation, gene, available, model in cases:
        original = a or b
        key = tuple(sorted([original['gene_a'], original['gene_b']]))
        for item in [a, b]:
            if item:
                validate_row(item)
        result = comparison(key, a, b)
        assert [result[f] for f in ['relation', 'same_lexical_gene', 'both_lexical_models_available', 'same_lexical_model']] == [relation, gene, available, model]
        output.append(result)
        if a:
            left.append(a)
        if b:
            b = b.copy()
            b['gene_a'], b['gene_b'] = b['gene_b'], b['gene_a']
            right.append(b)
    frame, matrices, totals = expected(pd.DataFrame(left), pd.DataFrame(right))
    pd.testing.assert_frame_equal(pd.DataFrame(output)[FIELDS], frame)
    assert totals['union_pairs'] == 9 and totals['shared_pairs'] == 7
    assert totals['sequence_eligible_in_both'] == 5 and totals['both_lexical_models_available'] == 3
    assert totals['same_lexical_model'] == 2  # Missing/missing is never model agreement.
    bad = row(2, ['A'], 'wrong')
    try:
        comparison(('T2_a', 'T2_b'), left[1], bad)
    except ValueError:
        pass
    else:
        raise AssertionError('Conflicting same-gene model accepted')
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        inventory = root / 'native'
        inventory.mkdir()
        for g, rows in [('profile', left), ('mafft', right)]:
            table(inventory / (g + '_sequence_first_references.tsv'), rows)
        sp, proof = root / 'native-plan.json', root / 'native-proof.json'
        write(sp, dict(output=str(inventory)))
        rp = inventory / 'receipt.json'
        write(rp, dict(status='complete_full_sequence_first_sister_reference_inventory_pending_readback',
                       plan_sha256=sha(sp), guides=[dict(guide=g, candidates=len(rows)) for g, rows in [('profile', left), ('mafft', right)]],
                       artifacts={p.name: sha(p) for p in inventory.iterdir()}))
        write(proof, dict(status='passed_full_sequence_first_sister_reference_native_readback', producer_receipt_sha256=sha(rp)))
        out = root / 'comparison'
        pp = root / 'plan.json'
        write(pp, dict(inventory=str(inventory), source_plan=str(sp), native_readback=str(proof), output=str(out),
                       pins={str(p): sha(p) for p in [*inventory.iterdir(), sp, proof]}))
        subprocess.run([sys.executable, 'scripts/compare_sequence_first_references.py', '--plan', str(pp)], check=True, capture_output=True)
        reader = [sys.executable, 'scripts/readback_sequence_first_reference_comparison.py', '--plan', str(pp)]
        subprocess.run([*reader, '--output', str(root / 'passed.json')], check=True, capture_output=True)
        exports = [out / 'sequence_first_guide_comparison.tsv', out / 'context_matrix.tsv', out / 'choice_matrix.tsv']
        originals = {p: p.read_bytes() for p in exports}
        receipt = out / 'receipt.json'
        native_receipt = receipt.read_bytes()
        for i, (path, field, value) in enumerate([
            (exports[0], 'same_lexical_model', '1'),  # Unmodeled/unmodeled agreement is forbidden.
            (exports[0], 'relation', 'disjoint_nearest_sequence_gene_sets'),
            (exports[0], 'profile_sequence_first_chosen_gene', 'wrong'),
            (exports[1], 'pairs', '999'), (exports[2], 'pairs', '999')]):
            for p, blob in originals.items():
                p.write_bytes(blob)
            receipt.write_bytes(native_receipt)
            frame = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
            frame.loc[0, field] = value
            frame.to_csv(path, sep='\t', index=False)
            r = json.loads(receipt.read_text())
            r['artifacts'][path.name] = sha(path)
            write(receipt, r)
            result = subprocess.run([*reader, '--output', str(root / f'false-{i}.json')], capture_output=True)
            assert result.returncode != 0 and not (root / f'false-{i}.json').exists()
    print('Passed9 known union/context/tie/missing/model/role cases, full synthetic handoff and5 rehashed false exports.')


if __name__ == '__main__':
    main()
