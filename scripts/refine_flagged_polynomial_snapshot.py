#!/usr/bin/env python3
"""Refine every fit in a frozen, independently verified review snapshot."""
import csv
import json
import subprocess
import time
from collections import Counter
from pathlib import Path
import numpy as np
import psutil
from refine_matched_ml_analytic import refine
from cached_matched_ml import profiled_ml
from matched_mixed_covariance import MatchedCovariance
from screen_duplication_alignment_reuse import sha


def main():
    pp = Path('metadata/flagged_polynomial_snapshot_refinement_plan_20260928.json')
    plan = json.loads(pp.read_text())
    bindings = {str(pp): sha(pp), **plan['pins']}
    launch = json.loads(Path(plan['audit_launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            p = psutil.Process(launch['pid'])
            if abs(p.create_time() - launch['created']) > .01 or p.status() == psutil.STATUS_ZOMBIE:
                break
            assert p.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(x.split('=', 1) for x in subprocess.check_output(['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    root = Path(plan['inputs'])
    receipt = json.loads((root / 'receipt.json').read_text())
    proof = json.loads(Path(plan['audit_proof']).read_text())
    assert proof['status'] == 'passed_full_flagged_polynomial_input_readback'
    assert proof['source_receipt_sha256'] == sha(root / 'receipt.json')
    bindings[plan['audit_proof']] = sha(plan['audit_proof'])
    bindings.update(receipt['source_hashes'])
    bindings.update({str(root / n): h for n, h in receipt['artifacts'].items()})
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    manifest = {r['recipe']['fit_input_id']: r for r in json.loads((root / 'input_manifest.json').read_text())}
    rows = list(csv.DictReader((root / 'flagged_tree_fits.tsv').open(), delimiter='\t'))
    assert len(rows) == 601 and len(manifest) == 449
    factors = {}
    for path in Path(plan['factors']).glob('*.npz'):
        with np.load(path) as arrays:
            factors[path.stem] = arrays['factor']
    assert len(factors) == 5
    out = Path(plan['output'])
    out.mkdir(exist_ok=False)
    counts = Counter()
    completed = []
    for index, row in enumerate(rows):
        identifier, tree = row['fit_input_id'], row['tree']
        original = json.loads(Path(row['source_path']).read_text())
        assert sha(row['source_path']) == row['source_sha256']
        assert original['fit_input_id'] == identifier and original['tree'] == tree
        baseline = original['payload']
        assert baseline['status'] == 'ml_candidate_requires_optimization_review'
        entry = manifest[identifier]
        with np.load(root / entry['path']) as arrays:
            numeric = arrays['matrix']
            background, family = arrays['background'], arrays['family']
            factor = factors[tree][arrays['pattern_rows']]
        active = np.ptp(numeric[:, 1:], axis=0) > 1e-12
        assert np.all(abs(numeric[:, 1:][:, ~active]) <= 1e-12)
        scales = np.std(numeric[:, 1:][:, active], axis=0)
        design = np.column_stack([np.ones(len(numeric)), numeric[:, 1:][:, active] / scales])
        response = numeric[:, 0]
        np.testing.assert_array_equal(active, baseline['active_covariates'])
        np.testing.assert_array_equal(scales, baseline['covariate_scales'])
        initial = np.asarray(baseline['log1p_ratios'])
        direct = profiled_ml(MatchedCovariance(background, family, factor, 1., *np.expm1(initial)), design, response)
        for name in ['negative_profiled_ml', 'beta', 'profiled_scale']:
            np.testing.assert_allclose(direct[name], baseline[name], rtol=1e-8, atol=1e-7)
        fitted = refine(background, family, factor, design, response, initial,
                        maximum_ratio=baseline['maximum_ratio'], maxiter=1000)
        conversion = np.r_[1., 1. / scales]
        fitted.update(raw_unit_beta=(np.asarray(fitted['beta']) * conversion).tolist(),
                      raw_unit_conditional_beta_covariance=(np.asarray(fitted['conditional_beta_covariance']) * conversion[:, None] * conversion[None, :]).tolist(),
                      covariate_scales=scales.tolist(), active_covariates=active.tolist())
        result = dict(fit_input_id=identifier, tree=tree, polynomial_degree=int(row['polynomial_degree']),
                      source_fit_path=row['source_path'], source_fit_sha256=row['source_sha256'],
                      input_sha256=entry['sha256'], payload=fitted)
        path = out / (identifier + '__' + tree + '.json')
        with path.open('x') as handle:
            json.dump(result, handle, indent=2, allow_nan=False)
            handle.write('\n')
        assert json.loads(path.read_text()) == result
        counts[fitted['status']] += 1
        completed.append(dict(fit_input_id=identifier, tree=tree, path=path.name, sha256=sha(path), status=fitted['status']))
        (out / 'state.json').write_text(json.dumps(dict(completed=len(completed), total=601, counts=dict(counts)), indent=2) + '\n')
        print('Refined snapshot fit', index + 1, '/601', fitted['status'], flush=True)
    verify()
    (out / 'fit_manifest.json').write_text(json.dumps(completed, indent=2) + '\n')
    result = dict(status='complete_flagged_polynomial_snapshot_refinements_pending_independent_readback',
                  fits=len(completed), counts=dict(counts), source_hashes=bindings, script_sha256=sha(__file__),
                  artifacts={p.name: sha(p) for p in out.iterdir()},
                  scope='Frozen 601-flag snapshot only; all original fits retained. Numerical refinements require independent output checks. Original full grid still needs completion/audit; no replacement, global optimum, calibrated inference or biological claim.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts']}), flush=True)


if __name__ == '__main__':
    main()
