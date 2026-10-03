"""Expose existing all-node draws in a separate future logger version.

The only model-program edit is a pure label view passed to the logger. No
distribution, density, transition, observation, initializer or draw changes.
"""
import hashlib
import math

from independent_short_sampler_outputs import NATIVE_ALPHABET, digest


BEFORE = ';let {catStates  = labeledNodeMap tree ancStates}'
AFTER = ';let {catStates  = labeledNodeMap (addInternalLabels tree) ancStates}'


def transform(source):
    assert source.count(BEFORE) == 1 and AFTER not in source
    assert source.count(';let {ancStates  = prop_anc_cat_states properties}') == 1
    assert source.count('ancestralAlignment tree alignment (getSMap smodel) aa ancStates') == 1
    assert source.count('observe sequenceData (phyloCTMC tree alignment smodel scale1)') == 1
    result = source.replace(BEFORE, AFTER)
    assert result.replace(AFTER, BEFORE) == source
    return result


def fresh_seeds(chain_ids, forbidden, namespace='fungal-full-node-logger-future-20261003-v1'):
    """Deterministic future namespace, disjoint from all supplied used seeds."""
    assert len(chain_ids) == len(set(chain_ids)) == 1620
    used = set(forbidden); seeds = {}
    for cid in sorted(chain_ids):
        counter = 0
        while True:
            raw = hashlib.sha256((namespace + ':' + cid + ':' + str(counter)).encode()).digest()
            seed = int.from_bytes(raw[:8], 'big') % (2**31 - 1) + 1
            if seed not in used: break
            counter += 1
        used.add(seed); seeds[cid] = seed
    assert len(set(seeds.values())) == 1620 and set(seeds.values()).isdisjoint(forbidden)
    return seeds


def validate_frame(frame, iteration, sequences, tips):
    """Exact all-node state/sequence correspondence, including unanchored residues."""
    assert set(frame) == {'iter', 'catStates', 'properties', 'conditions'}
    assert type(frame['iter']) is int and frame['iter'] == iteration
    assert set(frame['catStates']) == set(sequences)
    assert set(tips) < set(sequences) and frame['conditions'] == {}
    assert set(frame['properties']) == {'rate'}
    rates = frame['properties']['rate']; assert len(rates) == 4 and all(len(r) == 20 for r in rates)
    for row in rates:
        assert all(type(v) in [float,int] and math.isfinite(v) and v >= 0 for v in row)
        assert all(v == row[0] for v in row)
    assert math.isclose(sum(r[0] for r in rates)/4, 1, rel_tol=1e-10, abs_tol=1e-10)
    pairs = {'tip_pairs': 0, 'ancestral_pairs': 0}
    for name, sequence in sequences.items():
        value = frame['catStates'][name]; assert set(value) == {'categories', 'states'}
        cats, states = value['categories'], value['states']; letters = sequence.replace('-', '')
        assert len(cats) == len(states) == len(letters)
        assert all(type(v) is int and 0 <= v < 4 for v in cats)
        assert all(type(v) is int and 0 <= v < 20 for v in states)
        assert ''.join(NATIVE_ALPHABET[v] for v in states) == letters
        pairs['tip_pairs' if name in tips else 'ancestral_pairs'] += len(states)
    return dict(**pairs, native_nodes=len(sequences), native_tips=len(tips),
        native_ancestors=len(sequences)-len(tips), frame_sha256=digest(frame),
        ancestral_categories_available=True, posterior_qualified=False, scientific_eligibility=False)
