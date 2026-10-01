#!/usr/bin/env python3
"""Reject rehashed false exports against all actual crossed PMSF sources."""
import argparse
import contextlib
import csv
import io
import json
import shutil
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch
import readback_four_run_pmsf as reader
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); original_plan = json.loads(args.plan.read_text()); original_root = Path(original_plan['output'])
    original_r = json.loads((original_root / 'receipt.json').read_text()); rejected = []
    with tempfile.TemporaryDirectory(prefix='pmsf-full-false-exports-') as tmp:
        base = Path(tmp); root = base / 'export'; root.mkdir()
        for name in original_r['artifacts']: shutil.copyfile(original_root / name, root / name)
        plan = {**original_plan, 'output': str(root)}; pp = base / 'plan.json'; pp.write_text(json.dumps(plan, indent=2) + '\n')
        r = dict(original_r); r['plan_sha256'] = sha(pp); r['source_hashes'] = dict(original_r['source_hashes'])
        del r['source_hashes'][str(args.plan)]; r['source_hashes'][str(pp)] = sha(pp)
        rp = root / 'receipt.json'; rp.write_text(json.dumps(r, indent=2) + '\n')
        def audit(name):
            with patch.object(sys, 'argv', ['readback_four_run_pmsf.py', '--plan', str(pp), '--output', str(base / (name + '.json'))]), contextlib.redirect_stdout(io.StringIO()): reader.main()
        audit('baseline')
        originals = {name: (root / name).read_bytes() for name in original_r['artifacts']}
        changes = ['changed_RF', 'invented_consensus_SH_aLRT', 'promoted_weak_conflict', 'invalid_quartet',
                   'invented_role_boundary', 'missing_presence_cell', 'duplicate_comparison', 'missing_comparison', 'changed_branch_length']
        for name in changes:
            for file, content in originals.items(): (root / file).write_bytes(content)
            file = ('comparisons.tsv' if name in ['changed_RF', 'duplicate_comparison', 'missing_comparison'] else
                    'conflicts.tsv' if name in ['promoted_weak_conflict', 'invalid_quartet'] else
                    'rooting_boundary.tsv' if name == 'invented_role_boundary' else 'split_presence.tsv')
            with (root / file).open() as f:
                cr = csv.DictReader(f, delimiter='\t'); columns = cr.fieldnames; rows = list(cr)
            if name == 'changed_RF': rows[0]['rf_distance'] = str(int(rows[0]['rf_distance']) + 1)
            elif name == 'duplicate_comparison': rows.append(dict(rows[0]))
            elif name in ['missing_comparison', 'missing_presence_cell']: rows.pop()
            elif name == 'invented_consensus_SH_aLRT': next(r for r in rows if r['tree_type'] == 'consensus' and r['present'] == 'True')['sh_alrt'] = '99'
            elif name == 'promoted_weak_conflict': next(r for r in rows if r['both_support_criteria_met'] == 'False')['both_support_criteria_met'] = 'True'
            elif name == 'invalid_quartet': rows[0]['witness_quartet'] = ';'.join([rows[0]['witness_quartet'].split(';')[0]] * 4)
            elif name == 'invented_role_boundary': next(r for r in rows if r['boundary_split_present'] == 'False')['boundary_split_present'] = 'True'
            else: next(r for r in rows if r['present'] == 'True')['branch_length'] = '999'
            with (root / file).open('w') as f:
                cw = csv.DictWriter(f, fieldnames=columns, delimiter='\t', lineterminator='\n'); cw.writeheader(); cw.writerows(rows)
            r['artifacts'] = {file: sha(root / file) for file in original_r['artifacts']}; rp.write_text(json.dumps(r, indent=2) + '\n')
            try: audit(name)
            except (AssertionError, ValueError, KeyError): rejected.append(name)
            else: raise AssertionError('Accepted false full export: ' + name)
        assert len(rejected) == len(changes)
    result = dict(status='passed_full_four_run_pmsf_rehashed_false_export_checks', plan_sha256=sha(args.plan), producer_receipt_sha256=sha(original_root / 'receipt.json'),
                  actual_source_runs=4, actual_tree_views=8, actual_taxa=526, baseline_full_reader_passed=True, rejected_rehashed_false_exports=rejected,
                  checker_sha256=sha(__file__), scientific_eligibility=False,
                  scope='Complete actual data/schema/source gates exercised in isolated export copies. Source trees/audits and production outputs unchanged; alternative temporary-plan binding requalified before tests. Rehashed false table exports rejected independently. No synthetic taxa/source stubs or pilot, no scientific significance.')
    with args.output.open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
