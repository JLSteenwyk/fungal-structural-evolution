"""Prespecified mean-one inverse-reuse sensitivity controls, not error estimates."""
import numpy as np

SCHEMA = 'full-original-cohort-inverse-reuse-controls-v1'
POLICIES = ['uniform', 'background_node', 'background_pair', 'family_component']
ARRAYS = ['case_rows', 'reuse_counts', 'weights', 'reciprocal_diagonals']
CLAIMS = dict(calibrated_measurement_precision=False, inverse_probability_weighting=False,
              effective_sample_size_estimated=False, raw_reml_basis_qualification_complete=False,
              nonuniform_weighting_accepted=False, scientific_eligibility=False)


def calculate(keys):
    assert len(keys) == 3
    n = len(keys[0]); assert n > 0 and all(len(a) == n for a in keys)
    counts = np.empty((3, n), dtype=np.int64)
    weights = np.ones((4, n), dtype=np.float64)
    diagonals = weights.copy(); groups = []
    for j, values in enumerate(keys):
        assert np.asarray(values).ndim == 1
        labels, inverse, reuse = np.unique(values, return_inverse=True, return_counts=True)
        g = len(labels); groups.append(g); counts[j] = reuse[inverse]
        assert n * g < 2**53
        denominator = g * counts[j]
        weights[j + 1] = n / denominator
        diagonals[j + 1] = denominator / n
    assert np.isfinite(weights).all() and np.all(weights > 0)
    assert np.isfinite(diagonals).all() and np.all(diagonals > 0)
    return counts, weights, diagonals, groups


def policy_records(n, counts, weights, diagonals, groups):
    records = []
    for j, policy in enumerate(POLICIES):
        g = n if j == 0 else groups[j - 1]
        reuse = np.ones(n, dtype=np.int64) if j == 0 else counts[j - 1]
        uniform = bool(np.all(g * reuse == n))
        records.append(dict(policy=policy, represented_groups=g, records=n,
            maximum_reuse=int(reuse.max()), minimum_reuse=int(reuse.min()),
            weight_min=float(weights[j].min()), weight_max=float(weights[j].max()),
            diagonal_min=float(diagonals[j].min()), diagonal_max=float(diagonals[j].max()),
            weight_sum_float64=float(weights[j].sum()),
            exact_weight_sum=n, exact_group_total_numerator=n,
            exact_group_total_denominator=g, diagonal_is_exact_uniform_one=uniform,
            future_basis_disposition=('requires_uniform_named_fold_and_fresh_numerical_qualification'
                if uniform else 'requires_positive_diagonal_cone_and_fresh_numerical_qualification')))
    return records
