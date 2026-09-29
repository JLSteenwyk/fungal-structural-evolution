#!/usr/bin/env python3
"""Direct-likelihood finite differences for the first unresolved refined gradient."""
import csv
import json
from pathlib import Path
import numpy as np
from cached_matched_ml import profiled_ml, CachedMatchedML
from matched_mixed_covariance import MatchedCovariance
from matched_ml_gradient import evaluate_gradient
from screen_duplication_alignment_reuse import sha


def main():
    pp = Path('metadata/refinement_gradient_diagnostic_plan_20260928.json')
    plan = json.loads(pp.read_text())
    bindings = {str(pp): sha(pp), **plan['pins']}
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    source = json.loads(Path(plan['fit']).read_text())
    result = source['payload']
    assert result['status'] == 'ml_refinement_requires_review'
    with np.load(plan['inputs']) as a:
        matrix, bg, family, indices = a['matrix'], a['background'], a['family'], a['pattern_rows']
    with np.load(plan['factor']) as a:
        factor = a['factor'][indices]
    active = np.asarray(result['active_covariates'])
    scales = np.asarray(result['covariate_scales'])
    x = np.column_stack([np.ones(len(matrix)), matrix[:, 1:][:, active] / scales])
    y = matrix[:, 0]
    theta = np.asarray(result['log1p_ratios'])
    upper = np.log1p(result['maximum_ratio'])
    def objective(t):
        assert (t >= 0).all() and (t <= upper).all()
        return profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(t)), x, y)['negative_profiled_ml']
    base = objective(theta)
    np.testing.assert_allclose(base, result['negative_profiled_ml'], rtol=1e-10, atol=1e-8)
    analytic = evaluate_gradient(CachedMatchedML(bg, family, factor, x, y), np.expm1(theta))['log1p_ratio_gradient']
    np.testing.assert_allclose(analytic, result['analytic_gradient'], rtol=1e-10, atol=1e-10)
    rows = []
    for axis, name in enumerate(['background', 'family_component', 'species']):
        for step in [1e-4, 1e-5, 1e-6, 1e-7]:
            plus, minus = theta.copy(), theta.copy()
            plus[axis] += step
            minus[axis] -= step
            high, low = objective(plus), objective(minus)
            finite_difference = (high - low) / (2 * step)
            rows.append(dict(component=name, log1p_step=step, objective_plus=high, objective_minus=low,
                             analytic_gradient=float(analytic[axis]), central_difference=finite_difference,
                             difference_from_analytic=finite_difference - analytic[axis]))
    out = Path(plan['output'])
    out.mkdir(exist_ok=False)
    table = out / 'finite_differences.tsv'
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    actual = list(csv.DictReader(table.open(), delimiter='\t'))
    assert actual == [{k: str(v) for k, v in row.items()} for row in rows]
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    receipt = dict(status='complete_direct_likelihood_gradient_diagnostic', source_hashes=bindings,
                   script_sha256=sha(__file__), fit_input_id=source['fit_input_id'], tree=source['tree'],
                   records=len(y), candidate_theta=theta.tolist(), objective=base,
                   original_projected_gradient=result['projected_gradient'], original_checks=result['checks'],
                   evaluations=25, finite_difference_rows=len(rows), artifacts={table.name: sha(table)},
                   scope='One flagged fit, fixed parameter point, four finite-difference steps per component using direct likelihood. Diagnostic only: step-size truncation and floating-point cancellation remain visible. No optimization, status replacement or threshold relaxation.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(rows, indent=2), flush=True)


if __name__ == '__main__':
    main()
