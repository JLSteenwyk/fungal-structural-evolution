#!/usr/bin/env python3
"""Known sequence-first/availability choices plus end-to-end rehashed-export rejection."""
import csv
import json
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

from inventory_duplication_sister_references import tree_index, reference_summary
from inventory_sequence_first_sister_references import sequence_first
from readback_duplication_sister_references import index_tree, check_row
from readback_sequence_first_sister_references import reconstruct_sequence_first
from screen_duplication_alignment_reuse import sha


def save(path, record):
    Path(path).write_text(json.dumps(record) + '\n')


def write_rows(path, rows):
    with Path(path).open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    cases = [
        ('((F1_a:1,F1_b:2)n1:3,F2_r:4)n0;', {'F2_r': ('M', 6)}, (), 'sequence_first_lexical_reference_modeled'),
        ('((F1_a:1,F1_b:2)n1:3,(F2_r:1,F3_s:2)n2:3)n0;', {'F3_s': ('M', 10)}, (), 'nearest_sequence_ties_unmodeled_farther_reference_used'),
        ('((F1_a:1,F1_b:2)n1:3,(F2_r:1,F3_s:1)n2:3)n0;', {'F3_s': ('M', 10)}, (), 'lexical_reference_unmodeled_other_nearest_tie_modeled'),
        ('((F1_a:1,F1_b:2)n1:3,F2_r:4)n0;', {}, (), 'no_modeled_sister_reference'),
        ('((F1_a:1,F1_b:2)n1:3,F2_r:4)n0;', {'F2_r': ('M', 6)}, ('n0',), 'ineligible_parent_context'),
        ('((F1_a:1,F1_b:2)n1:3,F2_r:4,F3_s:4)n0;', {'F2_r': ('M', 6)}, (), 'ineligible_parent_context'),
        ('((F1_a:1,F1_b:2)n1:3,(F1_r:1,F3_s:2)n2:3)n0;', {'F3_s': ('M', 10)}, (), 'ineligible_parent_context'),
        ('(F1_a:1,F1_b:2)n1;', {}, (), 'ineligible_parent_context'),
        ('(((F1_a:1,F1_b:2)n1:3,F2_r:4)n0:1000000000000000,F3_s:1)n9;', {'F2_r': ('M', 6), 'F3_s': ('M', 10)}, (), 'sequence_first_lexical_reference_modeled'),
        ('((F1_a:1,F1_b:2)n1:3,(F2_r:1,F3_s:1)n2:3)n0;', {'F2_r': ('M', 6), 'F3_s': ('M', 10)}, (), 'sequence_first_lexical_reference_modeled'),
    ]
    results = []
    for i, (tree, models, duplications, status) in enumerate(cases):
        row = dict(family='OG' + str(i), gene_node='n1', taxon_id='F1', gene_a='F1_a', gene_b='F1_b')
        original = reference_summary(row, tree_index(tree), set(duplications), models, {'F1', 'F2', 'F3'})
        row.update(original)
        left = sequence_first(row, tree_index(tree), models)
        # Reader sees serialized native source values, just as in production.
        observed_row = {k: json.dumps(v) if isinstance(v, list) else str(v) for k, v in row.items()}
        right = reconstruct_sequence_first(observed_row, index_tree(tree), models)
        saved = {k: json.dumps(v) if isinstance(v, list) else str(v) for k, v in left.items()}
        assert check_row(saved, right) < 1e-12 and left['sequence_first_status'] == status
        results.append(left)
    assert results[1]['availability_distance_increment'] == 1.0
    assert results[2]['sequence_first_nearest_genes'] == ['F2_r', 'F3_s']
    assert results[2]['sequence_first_chosen_gene'] == 'F2_r' and results[2]['sequence_first_chosen_model'] == ''
    assert results[0]['sequence_first_minimum_distance'] == results[8]['sequence_first_minimum_distance'] == 7
    assert [r['version'] for r in results[9]['sequence_first_tied_model_assignments']] == [6, 10]

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / 'source'
        source.mkdir()
        bridge = root / 'bridge.sqlite'
        models = {'F3_s': ('M', 10)}
        with sqlite3.connect(bridge) as db:
            db.execute('CREATE TABLE structures(taxon_id,protein_id,model_id,version)')
            db.execute('INSERT INTO structures VALUES(?,?,?,?)', ('F3', 's', 'M', 10))
        trees = root / 'trees.txt'
        tree = cases[1][0]
        trees.write_text('OG: ' + tree + '\n')
        row = dict(family='OG', gene_node='n1', taxon_id='F1', gene_a='F1_a', gene_b='F1_b')
        row.update(reference_summary(row, tree_index(tree), set(), models, {'F1', 'F2', 'F3'}))
        for guide in ['profile', 'mafft']:
            write_rows(source / (guide + '_sister_references.tsv'), [row])
        review = root / 'review-plan.json'
        save(review, dict(bridge=str(bridge), guides=[dict(guide=g, trees=str(trees)) for g in ['profile', 'mafft']]))
        sp = root / 'source-plan.json'
        save(sp, dict(output=str(source), bridge=str(bridge), review_plan=str(review)))
        save(source / 'receipt.json', dict(status='complete_duplication_sister_reference_inventory', plan_sha256=sha(sp),
                                           guides=[dict(guide=g, candidates=1, counts={'provisional_reference_available': 1}) for g in ['profile', 'mafft']],
                                           artifacts={p.name: sha(p) for p in source.iterdir()}))
        native_proof = root / 'source-proof.json'
        save(native_proof, dict(status='passed_full_duplication_sister_reference_readback', producer_receipt_sha256=sha(source / 'receipt.json')))
        plan = root / 'plan.json'
        output = root / 'output'
        pins = {str(p): sha(p) for p in [sp, review, bridge, trees, native_proof, *source.iterdir()]}
        save(plan, dict(reference_plan=str(sp), references=str(source), reference_readback=str(native_proof), output=str(output), pins=pins))
        subprocess.run([sys.executable, 'scripts/inventory_sequence_first_sister_references.py', '--plan', str(plan)], check=True, capture_output=True)
        command = [sys.executable, 'scripts/readback_sequence_first_sister_references.py', '--plan', str(plan)]
        proof = root / 'proof.json'
        result = subprocess.run(command + ['--output', str(proof)], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        result = json.loads(proof.read_text())
        assert sum(r['candidates'] for r in result['guides']) == 2
        pristine = {p.name: p.read_bytes() for p in output.iterdir()}
        table = output / 'profile_sequence_first_references.tsv'

        def reject(number, field, value):
            rows = list(csv.DictReader(table.open(), delimiter='\t'))
            rows[0][field] = value
            write_rows(table, rows)
            r = json.loads((output / 'receipt.json').read_text())
            r['artifacts'][table.name] = sha(table)
            save(output / 'receipt.json', r)
            dest = root / ('rejected-' + str(number) + '.json')
            result = subprocess.run(command + ['--output', str(dest)], capture_output=True)
            assert result.returncode != 0 and not dest.exists()
            for name, data in pristine.items():
                (output / name).write_bytes(data)

        reject(1, 'sequence_first_chosen_gene', 'F3_s')
        reject(2, 'availability_distance_increment', '0')
        reject(3, 'sequence_first_tied_model_assignments', '[{"gene":"F2_r","model_id":"M","version":10}]')
        reject(4, 'sequence_first_status', 'sequence_first_lexical_reference_modeled')
    print('PASS: ten known choice/availability/tie/context/root cases, versions6/10, full two-guide handoff and four rehashed false export rejections.')


if __name__ == '__main__':
    main()
