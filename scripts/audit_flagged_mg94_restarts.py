#!/usr/bin/env python3
"""Read back the complete restart grid from saved logs and model declarations."""
import argparse
import json
import math
from pathlib import Path
from audit_genus_branch_parameter_profiles import declarations, sha, table, unchanged_model_text


def close(a, b):
    if not math.isfinite(float(a)) or not math.isfinite(float(b)) or abs(float(a)-float(b)) > 1e-6:
        raise ValueError(f'Numerical mismatch: {a}, {b}')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--restarts', type=Path, required=True)
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    receipt = json.loads((a.restarts/'receipt.json').read_text())
    plan = json.loads(a.plan.read_text())
    assert receipt['status'] == 'complete_flagged_unconstrained_multistart_diagnostics'
    assert sha(a.plan) == receipt['plan_sha256']
    assert sha(Path(__file__).with_name('restart_flagged_genus_profile_optima.py')) == receipt['script_sha256']
    assert sha(Path(__file__).with_name('audit_genus_branch_parameter_profiles.py')) == receipt['audit_helper_sha256']
    checked = 0
    for artifacts in [receipt['artifacts']] + [p['artifacts'] for p in receipt['proofs']]:
        for name, digest in artifacts.items():
            assert sha(a.restarts/name) == digest, name
            checked += 1
    rows = table(a.restarts/'restart_points.tsv')
    summaries = table(a.restarts/'case_summary.tsv')
    labels = {'unconstrained', 'factor_0p1', 'factor_0p25', 'factor_0p5', 'factor_1', 'factor_2', 'factor_4', 'factor_10'}
    expected = {(case, label) for case in plan['case_ids'] for label in labels}
    assert len(rows) == len(expected) == plan['fits'] == receipt['fits']
    assert {(r['case_id'], r['start_solution']) for r in rows} == expected
    assert len(receipt['proofs']) == len(expected)
    assert {(p['case_id'], p['start']) for p in receipt['proofs']} == expected
    max_error = 0
    for row in rows:
        work = a.restarts/row['case_id']/row['start_solution']
        lines = (work/'optimize.log').read_text().splitlines()
        opt = [x.split('\t')[1:] for x in lines if x.startswith('OPTIMUM\t')]
        assert len(opt) == 1
        close(opt[0][0], row['log_likelihood']); close(opt[0][1], row['log_likelihood'])
        pairs = [x.split('\t')[1:] for x in lines if x.startswith('PARAM\t')]
        fitted = {n: float(v) for n, v in pairs}
        assert len(fitted) == len(pairs)
        original = (work/'start.bf').read_text(); saved = (work/'fit.bf').read_text()
        values, fixed = declarations(saved); starts, start_fixed = declarations(original)
        assert not fixed and not start_fixed and values == fitted and set(starts) == set(values)
        assert all(math.isfinite(v) and v >= 0 for v in values.values())
        assert unchanged_model_text(saved, values) == unchanged_model_text(original, starts)
        omega = [v for n, v in values.items() if n.endswith('.omega')]
        assert len(omega) == 1; close(omega[0], row['omega'])
        replay = [float(x.split('=', 1)[1]) for x in (work/'readback.log').read_text().splitlines() if x.startswith('READBACK=')]
        assert len(replay) == 1
        error = replay[0]-float(row['log_likelihood'])
        close(error, row['fresh_readback_error']); close(error, 0)
        close(float(row['log_likelihood'])-float(row['starting_log_likelihood']), row['gain_over_start'])
        assert float(row['gain_over_start']) >= -1e-5
        max_error = max(max_error, abs(error))
    assert len(summaries) == len(plan['case_ids'])
    assert {r['case_id'] for r in summaries} == set(plan['case_ids'])
    for summary in summaries:
        selected = [r for r in rows if r['case_id'] == summary['case_id']]
        best = max(selected, key=lambda r: float(r['log_likelihood']))
        assert int(summary['starts']) == len(selected) and summary['best_start'] == best['start_solution']
        for field, source in [('best_log_likelihood', 'log_likelihood'), ('best_target_t', 'target_t'), ('best_omega', 'omega'), ('gain_over_previous_unconstrained', 'gain_over_previous_unconstrained'), ('gain_over_previous_best_evaluated', 'gain_over_previous_best_grid_or_unconstrained')]:
            close(summary[field], best[source])
        close(summary['range_of_restarted_log_likelihoods'], max(float(r['log_likelihood']) for r in selected)-min(float(r['log_likelihood']) for r in selected))
    result = {'status': 'passed_complete_restart_artifact_and_log_readback', 'cases': len(summaries), 'fits': len(rows), 'artifact_hashes_checked': checked, 'maximum_saved_readback_error': max_error, 'source_receipt_sha256': sha(a.restarts/'receipt.json'), 'script_sha256': sha(Path(__file__)), 'interpretation': 'Saved logs, parameter exports, unchanged model definitions and summary arithmetic agree. Does not rerun likelihoods, prove global optimality, validate biological assumptions or establish selection eligibility.'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
