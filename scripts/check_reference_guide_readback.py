#!/usr/bin/env python3
"""Check all guide-reference relations and reject a rehashed false export."""
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

from readback_duplication_reference_guide_comparison import compare_sources


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value) + '\n')


def row(n, refs, status='provisional_reference_available', version='6'):
    return dict(gene_a='T0_a' + str(n), gene_b='T0_b' + str(n), family='OG', gene_node='n' + str(n),
                status=status, chosen_reference_gene=refs[0] if refs else '',
                reference_model='M' if refs else '', reference_version=version if refs else '',
                nearest_reference_genes=json.dumps(refs))


def main():
    a = pd.DataFrame([row(0, ['T1_r', 'T2_r']), row(1, ['T1_r']), row(2, ['T1_r', 'T2_r']),
                      row(3, [], 'no_modeled_nonfocal_sister'), row(4, ['T1_r'])])
    b = pd.DataFrame([row(0, ['T1_r', 'T2_r']), row(1, ['T3_r']), row(2, ['T2_r', 'T3_r']),
                      row(3, ['T1_r']), row(5, ['T2_r'])])
    result, _, summary = compare_sources(a, b)
    assert result['relation'].tolist() == ['same_nearest_reference_set', 'disjoint_nearest_reference_sets',
                                         'overlapping_nearest_reference_sets', 'not_provisional_in_both',
                                         'profile_only_candidate', 'mafft_only_candidate']
    assert result['same_chosen_reference'].tolist() == ['1', '0', '0', '', '', '']
    assert summary['shared_duplicate_pairs'] == 4 and summary['provisional_in_both'] == 3 and len(result) == 6
    bad = b.copy()
    bad.loc[0, 'reference_version'] = '10'
    try:
        compare_sources(a, bad)
    except ValueError:
        pass
    else:
        raise AssertionError('Inconsistent frozen reference version accepted')
    try:
        compare_sources(pd.concat([a, a.iloc[[0]]]), b)
    except ValueError:
        pass
    else:
        raise AssertionError('Repeated source pair accepted')
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        inventory, review = root / 'inventory', root / 'review'
        inventory.mkdir()
        review.mkdir()
        for guide, frame in [('profile', a), ('mafft', b)]:
            frame.to_csv(inventory / (guide + '_sister_references.tsv'), sep='\t', index=False)
            frame[['gene_a', 'gene_b', 'family', 'gene_node']].to_csv(
                review / (guide + '_candidate_tree_checks.tsv'), sep='\t', index=False)
        save(review / 'receipt.json', dict(overlap={'shared_pairs': 4}, artifacts={p.name: sha(p) for p in review.iterdir()}))
        plan_path = root / 'plan.json'
        save(plan_path, dict(review=str(review), pins={str(review / 'receipt.json'): sha(review / 'receipt.json')}))
        guides = [dict(guide=g, candidates=5) for g in ['profile', 'mafft']]
        save(inventory / 'receipt.json', dict(status='complete_duplication_sister_reference_inventory',
                                             plan_sha256=sha(plan_path), guides=guides,
                                             artifacts={p.name: sha(p) for p in inventory.iterdir()}))
        proof = root / 'native-proof.json'
        save(proof, dict(status='passed_full_duplication_sister_reference_readback',
                         producer_receipt_sha256=sha(inventory / 'receipt.json'), guides=guides))
        output = root / 'comparison'
        subprocess.run([sys.executable, 'scripts/compare_duplication_sister_references.py',
                        '--inventory', str(inventory), '--plan', str(plan_path), '--review', str(review),
                        '--output', str(output)], check=True, capture_output=True)
        command = [sys.executable, 'scripts/readback_duplication_reference_guide_comparison.py',
                   '--source', str(output), '--inventory', str(inventory), '--reference-readback', str(proof),
                   '--review', str(review)]
        checked = root / 'checked.json'
        subprocess.run(command + ['--output', str(checked)], check=True, capture_output=True)
        assert json.loads(checked.read_text())['union_pairs'] == 6
        table_path = output / 'pair_reference_comparison.tsv'
        records = pd.read_csv(table_path, sep='\t', dtype=str, keep_default_na=False)
        records.loc[0, 'relation'] = 'disjoint_nearest_reference_sets'
        records.to_csv(table_path, sep='\t', index=False)
        receipt = json.loads((output / 'receipt.json').read_text())
        receipt['artifacts'][table_path.name] = sha(table_path)
        save(output / 'receipt.json', receipt)
        rejected = root / 'rejected.json'
        assert subprocess.run(command + ['--output', str(rejected)], capture_output=True).returncode != 0
        assert not rejected.exists()
    print('All six reference relations, complete outer join, ineligible blanks, frozen version/duplicate rejection, full saved-row replay and rehashed false relation rejection passed.')


if __name__ == '__main__':
    main()
