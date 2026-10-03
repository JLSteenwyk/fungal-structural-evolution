#!/usr/bin/env python3
"""Full categorical contracts against the locked oracle; no posterior acceptance."""
import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
import time
import warnings

import arviz
import numpy as np

from ancestral_categorical_diagnostics import diagnose_states as oracle
from ancestral_state_patterns import group_patterns as original_patterns
from independent_ancestral_categorical_diagnostics import compare_states, diagnose_states, group_patterns


def fixtures():
    rng = np.random.default_rng(20261004)
    for draws in [19, 20, 21, 50, 75, 500, 750]:
        yield f'constant-{draws}', np.zeros((4, draws), dtype=np.uint8)
        yield f'chain-constants-{draws}', np.repeat(np.arange(4, dtype=np.uint8)[:, None], draws, axis=1)
        yield f'binary-iid-{draws}', rng.integers(0, 2, (4, draws), dtype=np.uint8)
        yield f'all-labels-iid-{draws}', rng.integers(0, 22, (4, draws), dtype=np.uint8)
        yield f'unknown-only-{draws}', np.full((4, draws), 20, dtype=np.uint8)
        yield f'unknown-gap-{draws}', rng.integers(20, 22, (4, draws), dtype=np.uint8)
        halves = np.tile(np.r_[np.zeros(draws // 2, dtype=np.uint8),
            np.ones(draws - draws // 2, dtype=np.uint8)], (4, 1))
        yield f'constant-halves-{draws}', halves
        yield f'alternating-{draws}', np.tile(np.arange(draws, dtype=np.uint8) % 2, (4, 1))
        values = np.zeros((4, draws), dtype=np.uint8); values[0, -1] = 20
        yield f'one-rare-observation-{draws}', values
        values = rng.integers(0, 2, (4, draws), dtype=np.uint8); values[0] = 0
        yield f'one-constant-chain-{draws}', values
        values = rng.integers(0, 2, (4, draws), dtype=np.uint8); values[0] = 20
        yield f'chain-frequency-mismatch-{draws}', values
        values = rng.integers(0, 22, (4, draws), dtype=np.uint8)
        for i in range(1, draws):
            keep = rng.random(4) < .95; values[keep, i] = values[keep, i - 1]
        yield f'sticky-all-labels-{draws}', values


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert arviz.__version__ == '0.22.0' and np.__version__ == '2.2.6'
    started = time.perf_counter(); alphabet = list('ACDEFGHIKLMNPQRSTVWYX-')
    maxima = {}; reviews = []; states = Counter(); names = []; permutation_tv_error = 0.
    rng = np.random.default_rng(431)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        for name, values in fixtures():
            original, independent = oracle(values, alphabet), diagnose_states(values, alphabet)
            errors, unresolved = compare_states(original, independent, values, alphabet)
            for metric, error in errors.items(): maxima[metric] = max(maxima.get(metric, 0.), error)
            reviews.extend(dict(fixture=name, **review) for review in unresolved)
            states[original['status']] += 1; names.append(name)
            # Arbitrary simultaneous label/code permutation must preserve every
            # indicator, metric and chain distance after restoring label order.
            permutation = rng.permutation(len(alphabet))
            inverse = np.argsort(permutation)
            labels = [alphabet[s] for s in permutation]
            changed = diagnose_states(inverse[values], labels)
            assert changed['status'] == independent['status']
            assert set(changed['unobserved_states']) == set(independent['unobserved_states'])
            for left, right in zip(changed['pairwise_empirical_total_variation'],
                                   independent['pairwise_empirical_total_variation']):
                assert left['chains'] == right['chains']
                error = abs(left['total_variation'] - right['total_variation'])
                assert error <= 5e-16  # Permuted summation order may round differently.
                permutation_tv_error = max(permutation_tv_error, error)
            assert all(changed['indicators'][label] == independent['indicators'][label] for label in alphabet)
    invalid = [np.zeros((3, 50), dtype=int), np.zeros((4, 0), dtype=int),
        np.zeros((4, 50), dtype=float), np.full((4, 50), -1), np.full((4, 50), 22)]
    for values in invalid:
        try: diagnose_states(values, alphabet)
        except ValueError: pass
        else: raise AssertionError('Invalid state array accepted')
    for labels in [[], ['A', 'A']]:
        try: diagnose_states(np.zeros((4, 50), dtype=int), labels)
        except ValueError: pass
        else: raise AssertionError('Invalid alphabet accepted')
    values = rng.integers(0, 2, (4, 750), dtype=np.uint8)
    original = oracle(values, alphabet); independent = diagnose_states(values, alphabet)
    for dtype in [np.int8, np.int32, np.int64, np.uint8, np.uint32, np.uint64]:
        converted = values.astype(dtype)
        compare_states(oracle(converted, alphabet), diagnose_states(converted, alphabet), converted, alphabet)
    mutations = ['count', 'frequency', 'metric', 'status', 'missing-state', 'tv', 'unobserved', 'extra-field']
    for mutation in mutations:
        bad = copy.deepcopy(independent)
        if mutation == 'count': bad['indicators']['A']['counts_per_chain'][0] += 1
        elif mutation == 'frequency': bad['indicators']['A']['frequency_per_chain'][0] += .01
        elif mutation == 'metric': bad['indicators']['A']['screen']['bulk_ess'] += 1
        elif mutation == 'status': bad['status'] = 'observed_state_indicator_screens_pass_only' if original['status'] != 'observed_state_indicator_screens_pass_only' else 'categorical_mixing_requires_review'
        elif mutation == 'missing-state': del bad['indicators']['-']
        elif mutation == 'tv': bad['pairwise_empirical_total_variation'][0]['total_variation'] += .01
        elif mutation == 'unobserved': bad['unobserved_states'] = []
        elif mutation == 'extra-field': bad['scientific_eligibility'] = True
        try: compare_states(original, bad, values, alphabet)
        except AssertionError: pass
        else: raise AssertionError(('Altered export accepted', mutation))
    pattern_cases = []
    for draws in [1, 20, 50, 75]:
        for shape in [(1,), (17,), (4, 13), (2, 3, 5)]:
            values = rng.integers(0, 22, (4, draws, *shape), dtype=np.uint8)
            flat = values.reshape(4, draws, -1)
            if flat.shape[-1] > 2:
                flat[:, :, -1] = flat[:, :, 0]; flat[:, :, 1] = flat[:, :, 0]
            a, b = original_patterns(values); c, d = group_patterns(values)
            assert np.array_equal(a, c) and np.array_equal(b, d)
            assert np.array_equal(np.bincount(b.ravel()), np.bincount(d.ravel()))
            pattern_cases.append(dict(draws=draws, coordinate_shape=list(shape), patterns=len(a)))
    for values in [np.zeros((4, 2), dtype=np.uint8), np.zeros((4, 2, 0), dtype=np.uint8), np.zeros((4, 2, 3))]:
        try: group_patterns(values)
        except ValueError: pass
        else: raise AssertionError('Invalid pattern array accepted')
    paths = [__file__, 'scripts/independent_ancestral_categorical_diagnostics.py',
        'scripts/independent_ancestral_scalar_diagnostics_v2.py', 'scripts/ancestral_categorical_diagnostics.py',
        'scripts/ancestral_state_patterns.py', 'scripts/ancestral_chain_diagnostics.py']
    result = dict(status='passed_independent_ancestral_categorical_contracts',
        fixtures=len(names), fixture_names=names, fixture_status_counts=dict(states),
        alphabet=alphabet, every_declared_state_compared=True, label_permutation_cases=len(names),
        supported_integer_dtypes_compared=['int8', 'int32', 'int64', 'uint8', 'uint32', 'uint64'],
        maximum_absolute_errors=maxima, singular_indicator_reviews=reviews,
        maximum_permuted_total_variation_roundoff=permutation_tv_error,
        invalid_state_inputs_rejected=len(invalid) + 2, altered_exports_rejected=mutations,
        exact_pattern_cases=pattern_cases, invalid_pattern_inputs_rejected=3,
        arviz_oracle_version=arviz.__version__, numpy_version=np.__version__,
        source_hashes={str(p):sha(p) for p in paths}, elapsed_seconds=time.perf_counter()-started,
        scientific_eligibility=False, production_patterns_replayed=False,
        scope='Software fixtures only: separate counts, all-label indicator metrics, absent states, '
              'six descriptive distances, code/label permutation and lossless first-appearance '
              'pattern maps. No native parser validation, production posterior or biological pilot.')
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result))


if __name__ == '__main__': main()
