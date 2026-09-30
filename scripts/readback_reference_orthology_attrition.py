#!/usr/bin/env python3
"""Rebuild the complete reference attrition summaries with independent SQL aggregation."""
import argparse
import csv
import json
import sqlite3
from pathlib import Path

from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    plan = json.loads(args.plan.read_text())
    root = Path(plan['output'])
    source = Path(plan['source'])
    receipt = json.loads((root / 'receipt.json').read_text())
    source_receipt = json.loads((source / 'receipt.json').read_text())
    closure = json.loads(Path(plan['completion']).read_text())
    assert closure['status'] == 'complete_verified_full_reference_native_orthology' and len(closure['services']) == 2
    assert receipt['status'] == 'complete_full_reference_native_assignment_attrition_summary'
    assert receipt['plan_sha256'] == sha(args.plan)
    contexts = source / 'contexts.jsonl'
    assert receipt['target_contexts'] == closure['summary']['target_contexts'] == plan['target_contexts']
    bindings = {**receipt['source_hashes'], str(root / 'receipt.json'): sha(root / 'receipt.json'),
                str(contexts): source_receipt['artifacts'][contexts.name],
                **{str(root / k): v for k, v in receipt['artifacts'].items()}}
    assert closure['source_hashes'][str(contexts)] == bindings[str(contexts)]
    def verify():
        for path, digest in bindings.items(): assert sha(path) == digest, path
    verify()
    with sqlite3.connect(':memory:') as db:
        db.execute('CREATE TABLE contexts (guide TEXT,design TEXT,row_number INT,status TEXT,parent INT,genes INT,model INT,p TEXT,m TEXT,ties INT,both_ties INT,PRIMARY KEY(guide,design,row_number))')
        pending = []
        with contexts.open() as handle:
            for line in handle:
                row = json.loads(line)
                for design in ['availability', 'sequence_first']:
                    refs = row['designs'][design]
                    profiles = [(r['native_coorthology']['profile'], r['native_coorthology']['mafft']) for r in refs]
                    lexical = refs[0] if refs else None
                    assert not lexical or lexical['lexical_choice']
                    pending.append((row['source_guide'], design, row['source_row_number'], row['source']['status'],
                                    int(row['parent_context_eligible']), int(bool(refs)), int(bool(lexical and lexical['reference_model'])),
                                    profiles[0][0] if profiles else 'not_queried', profiles[0][1] if profiles else 'not_queried',
                                    len(refs), sum(p == m == 'both' for p, m in profiles)))
                if len(pending) >= 3000:
                    db.executemany('INSERT INTO contexts VALUES (?,?,?,?,?,?,?,?,?,?,?)', pending)
                    pending = []
        db.executemany('INSERT INTO contexts VALUES (?,?,?,?,?,?,?,?,?,?,?)', pending)
        n = db.execute('SELECT COUNT(*) FROM contexts').fetchone()[0]
        assert n == 2 * plan['target_contexts'] == receipt['context_design_records']
        assert db.execute('SELECT COUNT(*) FROM (SELECT guide,row_number FROM contexts GROUP BY guide,row_number HAVING COUNT(*) != 2)').fetchone()[0] == 0
        policies = {
            'all_source_contexts': '1',
            'parent_context_eligible': 'parent=1',
            'reference_gene_present': 'genes=1',
            'lexical_reference_model_present': 'model=1',
            'parent_eligible_reference_gene_present': 'parent=1 AND genes=1',
            'parent_eligible_lexical_model_present': 'parent=1 AND model=1',
            'parent_eligible_lexical_native_own_coorthology': "parent=1 AND ((guide='profile' AND p='both') OR (guide='mafft' AND m='both'))",
            'parent_eligible_lexical_native_both_guide_coorthology': "parent=1 AND p='both' AND m='both'",
            'parent_eligible_lexical_model_native_own_coorthology': "parent=1 AND model=1 AND ((guide='profile' AND p='both') OR (guide='mafft' AND m='both'))",
            'parent_eligible_lexical_model_native_both_guide_coorthology': "parent=1 AND model=1 AND p='both' AND m='both'",
            'parent_eligible_all_ties_native_both_guide_coorthology': 'parent=1 AND ties>0 AND both_ties=ties',
            'parent_eligible_any_tie_native_both_guide_coorthology': 'parent=1 AND both_ties>0',
        }
        expected = []
        for policy, clause in policies.items():
            for guide, design, count in db.execute('SELECT guide,design,SUM(CASE WHEN ' + clause + ' THEN 1 ELSE 0 END) FROM contexts GROUP BY guide,design'):
                expected.append(dict(guide=guide, design=design, policy=policy, contexts=str(count)))
        expected.sort(key=lambda r: (r['guide'], r['design'], r['policy']))
        with (root / 'context_policy_counts.tsv').open() as handle:
            assert list(csv.DictReader(handle, delimiter='\t')) == expected
        matrix = []
        columns = ['guide', 'design', 'source_context_status', 'parent_context_eligible', 'reference_gene_present',
                   'lexical_model_present', 'profile_native_coorthology', 'mafft_native_coorthology', 'contexts']
        for row in db.execute('SELECT guide,design,status,parent,genes,model,p,m,COUNT(*) FROM contexts GROUP BY guide,design,status,parent,genes,model,p,m ORDER BY guide,design,status,parent,genes,model,p,m'):
            matrix.append(dict(zip(columns, map(str, row))))
        with (root / 'lexical_reference_context_matrix.tsv').open() as handle:
            assert list(csv.DictReader(handle, delimiter='\t')) == matrix
        assert len(expected) == receipt['policy_rows'] and len(matrix) == receipt['matrix_rows']
        for guide, count in closure['summary']['guide_contexts'].items():
            assert db.execute('SELECT COUNT(*) FROM contexts WHERE guide=?', (guide,)).fetchone()[0] == 2 * count
    verify()
    proof = dict(status='passed_full_reference_native_assignment_attrition_sql_readback',
                 plan_sha256=sha(args.plan), producer_receipt_sha256=sha(root / 'receipt.json'),
                 target_contexts=plan['target_contexts'], context_design_records=n, policy_rows=len(expected),
                 matrix_rows=len(matrix), source_hashes=bindings, checker_sha256=sha(__file__), scientific_eligibility=False,
                 scope='Every completed context/design loaded into a uniquely keyed SQL table; all parent/missing-model/gene '
                       'states, lexical native coorthology and any/all tie sensitivities independently aggregated with SQL '
                       'and compared field-for-field to both complete producer tables. No producer summary code shared. '
                       'Native membership/journal source closure required. Counts overlap across guides/designs/policies '
                       'and are not biological orthology validations, independent events or structural effects.')
    args.output.write_text(json.dumps(proof, indent=2) + '\n')
    print(json.dumps({k: v for k, v in proof.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
