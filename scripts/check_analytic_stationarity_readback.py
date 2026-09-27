"""Exercise boundary signs, exact threshold and error retention in the readback."""
from copy import deepcopy
import json
import math
from pathlib import Path
from readback_full_analytic_stationarity import check_row, sha


def main():
    payload = dict(status='candidate_requires_optimization_review',
        log1p_ratios=[0., 1., math.log1p(10000)], maximum_ratio=10000,
        checks=dict(best_optimizer_success=True, upper_bound_contact=True,
                    all_full_face_starts_agree=True, projected_gradient_pass=False))
    row = dict(original_status=payload['status'], assessment='analytic_score_evaluated',
        analytic_gradient=[4., .001, -3.], analytic_projected_gradient=[0., .001, 0.],
        analytic_threshold_pass=True, other_optimizer_checks_pass=False,
        original_gradient_pass=False, all_numerical_checks_with_analytic_gradient=False)
    assert check_row(row, payload) == ['upper_bound_contact']
    reverse = deepcopy(row)
    reverse.update(analytic_gradient=[-4., .001, 3.],
                   analytic_projected_gradient=[-4., .001, 3.], analytic_threshold_pass=False)
    assert check_row(reverse, payload) == ['analytic_gradient_threshold', 'upper_bound_contact']
    failed = deepcopy(payload)
    failed['checks'].update(best_optimizer_success=False, all_full_face_starts_agree=False)
    assert check_row(row, failed) == ['optimizer_termination', 'upper_bound_contact', 'full_face_start_disagreement']
    error_payload = dict(status='fit_error_requires_review', error_type='ArithmeticError')
    error_row = dict(original_status='fit_error_requires_review',
                     assessment='original_fit_error_unassessed', error_type='ArithmeticError')
    assert check_row(error_row, error_payload) == ['original_fit_error']
    unresolved = dict(original_status=payload['status'], assessment='analytic_diagnostic_unresolved',
                      error_type='ArithmeticError', error='cancellation')
    assert check_row(unresolved, payload) == ['analytic_diagnostic_unresolved']
    mutations = {
        'wrong_boundary_projection': {'analytic_projected_gradient': [4., .001, -3.]},
        'lost_upper_bound_flag': {'other_optimizer_checks_pass': True},
        'false_clearance': {'all_numerical_checks_with_analytic_gradient': True},
        'changed_original_flag': {'original_gradient_pass': True},
        'threshold_mismatch': {'analytic_threshold_pass': False},
        'nonfinite_derivative': {'analytic_gradient': [float('nan'), .001, -3.]},
        'missing_dimension': {'analytic_gradient': [4., .001]},
        'changed_original_status': {'original_status': 'candidate_passed_numerical_optimization_checks'},
        'unknown_assessment': {'assessment': 'passed'},
    }
    rejected = []
    for name, change in mutations.items():
        bad = dict(row, **change)
        try:
            check_row(bad, payload)
        except ValueError:
            rejected.append(name)
        else:
            raise AssertionError('Accepted corruption: ' + name)
    try:
        check_row(dict(error_row, assessment='analytic_score_evaluated'), error_payload)
    except ValueError:
        rejected.append('lost_original_fit_error')
    else:
        raise AssertionError('Lost original fit error')
    result = dict(status='passed_analytic_stationarity_readback_boundary_and_rejection_checks',
                  valid_cases=5, rejected_corruptions=rejected,
                  checker_sha256=sha('scripts/readback_full_analytic_stationarity.py'),
                  script_sha256=sha(__file__),
                  scope='Scalar row checks only; full-grid accounting remains pending the production and analytic assessment.')
    with Path('metadata/analytic_stationarity_readback_checks_20260927.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
