"""Separate direct-lag implementation of the locked scalar-screen estimators.

Rank/folded split R-hat and bulk/tail ESS follow Vehtari et al. (2021),
https://arxiv.org/abs/1903.08008. Finite-lag ESS refinement and mean MCSE target
the declared ArviZ 0.22.0 convention. No ArviZ or project diagnostic imports.
These marginal estimators do not qualify a joint ancestral posterior.
"""
import numpy as np
from scipy.special import ndtri
from scipy.stats import rankdata


def split(values):
    half = values.shape[1] // 2
    return np.concatenate((values[:, :half], values[:, -half:]), axis=0)


def normal_scores(values):
    ranks = rankdata(values.ravel(), method='average').reshape(values.shape)
    return ndtri((ranks - 0.375) / (values.size + 0.25))


def variance_parts(values):
    draws = values.shape[1]
    within = np.var(values, axis=1, ddof=1).mean()
    between_mean = np.var(values.mean(axis=1), ddof=1)
    return within, within * (draws - 1) / draws + between_mean


def potential_scale(values):
    within, pooled = variance_parts(values)
    if within == 0:
        return float('inf') if pooled > 0 else float('nan')
    return float(np.sqrt(pooled / within))


def effective_size(values):
    """Direct dot-product lags, paired positive prefix and cumulative minimum.

    Unlike the FFT-based oracle this computes autocovariances directly. The
    last evaluated even lag is retained with the locked finite-lag refinement;
    the lower autocorrelation-time bound permits antithetical ESS > draw count.
    """
    chains, draws = values.shape
    total = chains * draws
    if np.ptp(values) < np.finfo(float).resolution:
        return float(total)
    centered = values - values.mean(axis=1, keepdims=True)
    acov = np.asarray([np.correlate(row, row, mode='full')[draws - 1:] / draws
                       for row in centered])
    within = acov[:, 0].mean() * draws / (draws - 1)
    pooled = within * (draws - 1) / draws + np.var(values.mean(axis=1), ddof=1)
    if not np.isfinite(pooled) or pooled <= 0:
        return float('nan')
    rho = 1 - (within - acov.mean(axis=0)) / pooled
    rho[0] = 1.
    # Pairs are (rho[0]+rho[1]), (rho[2]+rho[3]), ... . The
    # penultimate lag limit reproduces the explicitly declared finite horizon.
    limit = (draws - 3) // 2
    pairs = rho[:2 * (limit + 1)].reshape(-1, 2).sum(axis=1)
    if pairs[0] <= 0:
        last = 0
    else:
        stopped = np.flatnonzero(pairs[1:] <= 0)
        last = min(int(stopped[0]) + 1, limit) if len(stopped) else limit
    accepted = np.minimum.accumulate(pairs[:last])
    tail_even = rho[2 * last] if pairs[last] >= 0 or rho[2 * last] > 0 else 0.
    time = -1 + 2 * accepted.sum() + tail_even
    time = max(float(time), 1 / np.log10(total))
    return float(total / time) if np.isfinite(rho).all() else float('nan')


def diagnose(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[0] != 4:
        raise ValueError('Exactly four chains required')
    draws = values.shape[1]
    if draws < 20:
        return dict(status='insufficient_retained_draws', draws_per_chain=draws)
    if not np.isfinite(values).all():
        return dict(status='nonfinite_draws', draws_per_chain=draws)
    if any(np.all(row == row[0]) for row in values):
        return dict(status='constant_chain_requires_review', draws_per_chain=draws)
    halves = split(values)
    normalized = normal_scores(halves)
    folded = normal_scores(np.abs(halves - np.median(halves)))
    # The locked oracle uses ordered max: a globally constant folded score
    # contributes no defined scale check, while the bulk check remains used.
    rhat = float(max(potential_scale(normalized), potential_scale(folded)))
    bulk = effective_size(normalized)
    tails = [effective_size(split((values <= np.quantile(values, p)).astype(float))) for p in [.05, .95]]
    tail = float(np.min(tails))
    mean_ess = effective_size(halves)
    mcse = float(np.std(values, ddof=1) / np.sqrt(mean_ess))
    metrics = dict(rhat=rhat, bulk_ess=bulk, tail_ess=tail, mean_mcse=mcse)
    finite = all(np.isfinite(x) for x in metrics.values())
    passed = finite and rhat < 1.01 and min(bulk, tail) >= 400
    return dict(status='passes_scalar_screen_only' if passed else 'scalar_mixing_requires_review',
        draws_per_chain=draws, **{k: v if np.isfinite(v) else None for k, v in metrics.items()})


def compare(original, independent, values=None, atol=1e-8, rtol=1e-8):
    assert original.keys() == independent.keys()
    assert original['status'] == independent['status']
    assert original['draws_per_chain'] == independent['draws_per_chain']
    errors = {}; unresolved = []
    for name in ['rhat', 'bulk_ess', 'tail_ess', 'mean_mcse']:
        if name not in original:
            continue
        a, b = original[name], independent[name]
        if name == 'rhat' and values is not None:
            halves = split(np.asarray(values, dtype=float))
            if all(np.all(row == row[0]) for row in halves):
                assert original['status'] == independent['status'] == 'scalar_mixing_requires_review'
                assert (a is None or a > 1.01) and (b is None or b > 1.01)
                unresolved.append('rhat_zero_split_within_variance_requires_review')
                continue
        if (a is None) != (b is None):
            # Exactly constant split trajectories have zero mathematical
            # within-chain variance. Floating summation can turn the oracle's
            # undefined/infinite R-hat into a huge finite number. Preserve
            # this as unresolved numerical review, never a matching value.
            raise AssertionError(('Unexplained finite/undefined metric mismatch', name, a, b))
        if a is None:
            continue
        assert np.isfinite(a) and np.isfinite(b)
        assert abs(a - b) <= atol + rtol * abs(a), (name, a, b)
        errors[name] = abs(a - b)
    return errors, unresolved
