"""Exact integer/Fraction certificate for binary Geyer truncation boundaries.

This does not turn a mismatch into numerical agreement. It enumerates possible
roundoff decisions only at an exactly zero reachable autocorrelation pair.
"""
from fractions import Fraction
import math

import numpy as np

from independent_ancestral_scalar_diagnostics_v2 import split, type7_quantile, compare


def candidates(values):
    values = np.asarray(values)
    assert values.ndim == 2 and values.shape[0] >= 2 and values.shape[1] >= 4
    assert np.all((values == 0) | (values == 1))
    values = values.astype(np.int64); chains, draws = values.shape; total = chains * draws
    counts = [int(row.sum()) for row in values]
    if len(set(values.ravel().tolist())) == 1:
        return dict(ess=[float(total)], exact_zero_pair_indices=[])
    def acov(lag):
        numerator = 0
        for row,count in zip(values,counts):
            left, right = row[:draws-lag], row[lag:]
            numerator += (draws**2 * int(np.dot(left,right))
                - draws * count * int(left.sum() + right.sum()) + (draws-lag) * count**2)
        return Fraction(numerator, chains * draws**3)
    zero = acov(0); within = zero * draws / (draws - 1)
    between = Fraction(chains * sum(c*c for c in counts) - sum(counts)**2,
                       chains * (chains - 1) * draws**2)
    pooled = zero + between; assert pooled > 0
    rho = [Fraction(1)] + [1 - (within - acov(lag)) / pooled for lag in range(1,draws)]
    limit = (draws - 3) // 2
    pairs = [rho[2*i] + rho[2*i+1] for i in range(limit+1)]
    terminals = []; zeros = []
    for i,pair in enumerate(pairs):
        if pair == 0:
            zeros.append(i); terminals.append((i,True))
            if rho[2*i] < 0: terminals.append((i,False))
        elif pair < 0:
            terminals.append((i,False)); break
        if i == limit: terminals.append((i, pair >= 0))
    outcomes = []
    for last,stored in terminals:
        prefix = Fraction(0); running = pairs[0]
        for pair in pairs[:last]:
            running = min(running,pair); prefix += running
        even = rho[2*last] if stored or rho[2*last] > 0 else Fraction(0)
        tau = max(float(-1 + 2*prefix + even), 1 / math.log10(total))
        ess = total / tau
        if ess not in outcomes: outcomes.append(ess)
    assert outcomes
    return dict(ess=outcomes, exact_zero_pair_indices=zeros, split_chains=chains,
        split_draws=draws, counts=counts,
        scope='Only exact reachable zero pairs permit alternate stop/storage decisions; genuine negative pairs stop enumeration.')


def compare_with_binary_certificate(original, independent, values, atol=1e-8, rtol=1e-8):
    """Keep all ordinary checks strict; mark certified differing metrics unresolved."""
    try: return compare(original,independent,values=values,atol=atol,rtol=rtol)
    except AssertionError: pass
    values = np.asarray(values,dtype=float); levels = np.unique(values)
    assert len(levels) == 2 and np.isfinite(levels).all()
    binary = (values == levels[1]).astype(np.int64)
    bulk = candidates(split(binary)); mean = bulk
    tails = [candidates(split((values <= type7_quantile(values,p)).astype(np.int64))) for p in [.05,.95]]
    options = dict(bulk_ess=bulk['ess'], tail_ess=sorted(set(min(a,b) for a in tails[0]['ess'] for b in tails[1]['ess'])),
        mean_mcse=[float(np.std(values,ddof=1)/np.sqrt(ess)) for ess in mean['ess']])
    certificates = dict(bulk_ess=[bulk],tail_ess=tails,mean_mcse=[mean])
    left, right = dict(original), dict(independent); unresolved = []
    for metric in options:
        a,b = left.get(metric),right.get(metric)
        if a is None or b is None or abs(a-b) <= atol + rtol*abs(a): continue
        assert any(c['exact_zero_pair_indices'] for c in certificates[metric])
        assert all(any(abs(value-choice) <= atol + rtol*abs(value) for choice in options[metric]) for value in [a,b])
        unresolved.append(dict(metric=metric,reason='exact_binary_zero_pair_truncation_requires_review',
            original=a,independent=b,admissible_roundoff_outcomes=options[metric],certificates=certificates[metric]))
        # Exclude this metric from the ordinary-agreement error map, keeping
        # all source/disposition/definedness/Rhat and other metrics unchanged.
        left[metric] = right[metric] = None
    assert unresolved, 'No certified numerical discrepancy'
    errors, existing_reviews = compare(left,right,values=values,atol=atol,rtol=rtol)
    return errors, [*existing_reviews,*unresolved]
