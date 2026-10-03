"""Separate unordered-state accounting and direct-lag marginal diagnostics.

Imports no original categorical, pattern or ArviZ implementation. This is a
software verification backend, not qualification of an ancestral posterior.
"""
from itertools import combinations

import numpy as np

from independent_ancestral_scalar_diagnostics_v2 import diagnose
from binary_ess_boundary_certificate import compare_with_binary_certificate as compare


def diagnose_states(values, alphabet):
    values = np.asarray(values)
    if values.ndim != 2 or values.shape[0] != 4 or values.shape[1] == 0:
        raise ValueError('Four nonempty equal-length chains required')
    if not alphabet or len(set(alphabet)) != len(alphabet):
        raise ValueError('Distinct declared labels required')
    if (not np.issubdtype(values.dtype, np.integer)
            or np.any(values < 0) or np.any(values >= len(alphabet))):
        raise ValueError('Integer codes within the declared alphabet required')
    draws = values.shape[1]
    counts = np.stack([np.bincount(row.astype(np.intp), minlength=len(alphabet)) for row in values])
    frequencies = counts / draws
    indicators = {}
    for state, label in enumerate(alphabet):
        total = int(counts[:, state].sum())
        if total == 0:
            screen = dict(status='state_not_observed_probability_unresolved', draws_per_chain=draws)
        elif total == values.size:
            screen = dict(status='state_constant_in_all_chains_requires_review', draws_per_chain=draws)
        else:
            screen = diagnose(np.equal(values, state).astype(float))
        indicators[label] = dict(counts_per_chain=counts[:, state].tolist(),
            frequency_per_chain=frequencies[:, state].tolist(), screen=screen)
    observed = np.flatnonzero(counts.sum(axis=0))
    if len(observed) == 1:
        status = 'no_observed_state_variation_requires_review'
    elif all(indicators[alphabet[s]]['screen']['status'] == 'passes_scalar_screen_only'
             for s in observed):
        status = 'observed_state_indicator_screens_pass_only'
    else:
        status = 'categorical_mixing_requires_review'
    distances = [dict(chains=[a, b], total_variation=float(
        np.abs(frequencies[a] - frequencies[b]).sum() / 2))
        for a, b in combinations(range(4), 2)]
    return dict(status=status, draws_per_chain=draws, observed_states=len(observed),
        unobserved_states=[alphabet[s] for s in range(len(alphabet)) if s not in observed],
        indicators=indicators, pairwise_empirical_total_variation=distances,
        scope='Separate unordered-label counts and direct-lag indicator screens; '
              'unseen probabilities, joint mixing and model qualification remain unresolved.')


def compare_states(original, independent, values, alphabet, atol=1e-8, rtol=1e-8):
    """Compare every state, disposition and descriptive distance, not just totals."""
    fields = {'status', 'draws_per_chain', 'observed_states', 'unobserved_states',
              'indicators', 'pairwise_empirical_total_variation', 'scope'}
    assert set(original) == set(independent) == fields
    for name in ['status', 'draws_per_chain', 'observed_states', 'unobserved_states']:
        assert original[name] == independent[name], name
    assert list(original['indicators']) == list(independent['indicators']) == list(alphabet)
    maximum = {}; unresolved = []
    for state, label in enumerate(alphabet):
        left, right = original['indicators'][label], independent['indicators'][label]
        assert set(left) == set(right) == {'counts_per_chain', 'frequency_per_chain', 'screen'}
        for name in ['counts_per_chain', 'frequency_per_chain']:
            assert left[name] == right[name], (label, name)
        errors, reviews = compare(left['screen'], right['screen'],
            values=np.equal(values, state).astype(float), atol=atol, rtol=rtol)
        for name, error in errors.items():
            maximum[name] = max(maximum.get(name, 0.), error)
        unresolved.extend(dict(state=label, **review) if isinstance(review, dict) else
            dict(state=label, reason=review) for review in reviews)
    assert original['pairwise_empirical_total_variation'] == independent['pairwise_empirical_total_variation']
    return maximum, unresolved


def group_patterns(values):
    """Sort full temporal byte rows, then restore first-coordinate appearance order.

    Unlike the original per-coordinate hash dictionary, this uses lexicographic
    unique rows. Every chain/draw/state remains intact; labels are never ranked.
    """
    values = np.asarray(values)
    if values.dtype != np.uint8 or values.ndim < 3 or any(n == 0 for n in values.shape):
        raise ValueError('Nonempty uint8 chain/draw/coordinate array required')
    chains, draws = values.shape[:2]
    coordinates = int(np.prod(values.shape[2:]))
    if coordinates > np.iinfo(np.int32).max:
        raise ValueError('Coordinate map exceeds int32 capacity')
    rows = values.reshape(chains, draws, coordinates).transpose(2, 0, 1).reshape(coordinates, -1)
    _, first, inverse = np.unique(rows, axis=0, return_index=True, return_inverse=True)
    order = np.argsort(first, kind='stable')
    renumber = np.empty(len(first), dtype=np.int32)
    renumber[order] = np.arange(len(first), dtype=np.int32)
    patterns = rows[first[order]].reshape(len(first), chains, draws).copy()
    mapping = renumber[inverse].reshape(values.shape[2:])
    reconstructed = patterns[mapping.ravel()].transpose(1, 2, 0).reshape(values.shape)
    assert np.array_equal(reconstructed, values)
    return patterns, mapping
