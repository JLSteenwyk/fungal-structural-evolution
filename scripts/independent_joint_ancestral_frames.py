"""Decode same-record ancestral sequence/category draws and residue coordinates.

No joint logger validation or original category projection is imported. This
reader checks every node before exporting candidate projections. It does not
qualify posterior mixing, native likelihoods or full-grid execution.
"""
import math
from pathlib import Path

import numpy as np

from independent_native_ancestral_alignment import fasta_records
from independent_short_sampler_outputs_v2 import NATIVE_ALPHABET, array_digest, digest


GAP = 255
ARRAY_NAMES = {'states', 'categories', 'anchor_columns', 'anchor_tip_indices',
               'anchor_tip_offsets', 'unanchored_candidate_indices',
               'unanchored_columns', 'unanchored_states', 'unanchored_categories',
               'category_rates'}


def decode(frame, iteration, observed, runtime_labels, candidates):
    """Candidate rows sorted by source node; anchors sorted by tip then residue.

    Anchor columns can repeat: different tips can anchor the same native column.
    Unanchored residues are exported once per candidate/native column. Gap 255
    is outside both native supports; no category is assigned to a gap.
    """
    assert set(frame) == {'iter', 'catStates', 'alignmentLines', 'properties', 'conditions'}
    assert type(frame['iter']) is int and frame['iter'] == iteration
    assert type(iteration) is int and iteration >= 0
    assert observed and set(observed) < set(runtime_labels)
    assert candidates and len(set(candidates.values())) == len(candidates)
    assert set(candidates.values()) <= set(runtime_labels) - set(observed)
    lines = frame['alignmentLines']
    assert type(lines) is list and lines
    assert all(type(line) is str and line and not any(ord(c) < 32 or c in '\"\\' for c in line) for line in lines)
    # Require exact identifiers and uppercase sequence tokens before the shared
    # FASTA tokenizer can normalize whitespace/case.
    for line in lines:
        if line.startswith('>'):
            assert line[1:] and not any(c.isspace() for c in line[1:])
        else:
            assert all(c in NATIVE_ALPHABET + '-' for c in line)
    sequences = fasta_records(lines)
    assert set(sequences) == set(runtime_labels) == set(frame['catStates'])
    widths = {len(sequence) for sequence in sequences.values()}
    assert len(widths) == 1
    width = widths.pop()
    assert 0 < width <= np.iinfo(np.int32).max
    assert frame['conditions'] == {} and set(frame['properties']) == {'rate'}
    rates = frame['properties']['rate']
    assert type(rates) is list and len(rates) == 4
    for row in rates:
        assert type(row) is list and len(row) == 20
        assert all(type(v) in (int, float) and math.isfinite(v) and v >= 0 for v in row)
        assert all(v == row[0] for v in row)
    assert math.isclose(sum(row[0] for row in rates) / 4, 1, rel_tol=1e-10, abs_tol=1e-10)
    tip_pairs = ancestral_pairs = unknowns = 0
    for label, sequence in sequences.items():
        record = frame['catStates'][label]
        assert set(record) == {'states', 'categories'}
        states, cats = record['states'], record['categories']
        assert type(states) is list and type(cats) is list
        letters = sequence.replace('-', '')
        assert len(states) == len(cats) == len(letters)
        assert all(type(v) is int and 0 <= v < 20 for v in states)
        assert all(type(v) is int and 0 <= v < 4 for v in cats)
        assert ''.join(NATIVE_ALPHABET[v] for v in states) == letters
        if label in observed:
            source = observed[label]
            assert source and all(c in NATIVE_ALPHABET + 'X' for c in source)
            assert len(source) == len(letters)
            assert all(a == 'X' or a == b for a, b in zip(source, letters))
            unknowns += source.count('X'); tip_pairs += len(states)
        else:
            ancestral_pairs += len(states)
    tips = sorted(observed); nodes = sorted(candidates)
    assert len(tips) <= np.iinfo(np.int32).max and len(nodes) <= np.iinfo(np.int32).max
    anchor_columns = []; tip_indices = []; tip_offsets = []
    occupied = np.zeros(width, dtype=bool)
    for index, tip in enumerate(tips):
        columns = [col for col, c in enumerate(sequences[tip]) if c != '-']
        anchor_columns.extend(columns)
        tip_indices.extend([index] * len(columns)); tip_offsets.extend(range(len(columns)))
        occupied[columns] = True
    columns = np.asarray(anchor_columns, dtype=np.int32)
    states = np.full((len(nodes), len(columns)), GAP, dtype=np.uint8)
    categories = np.full_like(states, GAP)
    free_node = []; free_col = []; free_state = []; free_cat = []
    for index, node in enumerate(nodes):
        label = candidates[node]; record = frame['catStates'][label]
        full_states = np.full(width, GAP, dtype=np.uint8)
        full_cats = np.full(width, GAP, dtype=np.uint8)
        residue = 0
        for col, letter in enumerate(sequences[label]):
            if letter == '-': continue
            state = record['states'][residue]; cat = record['categories'][residue]
            full_states[col] = state; full_cats[col] = cat
            if not occupied[col]:
                free_node.append(index); free_col.append(col)
                free_state.append(state); free_cat.append(cat)
            residue += 1
        assert residue == len(record['states'])
        states[index] = full_states[columns]; categories[index] = full_cats[columns]
    arrays = dict(states=states, categories=categories, anchor_columns=columns,
        anchor_tip_indices=np.asarray(tip_indices, dtype=np.int32),
        anchor_tip_offsets=np.asarray(tip_offsets, dtype=np.int32),
        unanchored_candidate_indices=np.asarray(free_node, dtype=np.int32),
        unanchored_columns=np.asarray(free_col, dtype=np.int32),
        unanchored_states=np.asarray(free_state, dtype=np.uint8),
        unanchored_categories=np.asarray(free_cat, dtype=np.uint8),
        category_rates=np.asarray([row[0] for row in rates], dtype=np.float64))
    assert np.array_equal(states == GAP, categories == GAP)
    summary = dict(iteration=iteration, frame_sha256=digest(frame),
        alignment_sha256=digest(sequences), native_nodes=len(sequences),
        native_tips=len(tips), native_ancestors=len(sequences)-len(tips),
        tip_pairs=tip_pairs, ancestral_pairs=ancestral_pairs,
        unknown_observed_tip_positions=unknowns, candidate_source_nodes=nodes,
        candidate_runtime_nodes=[candidates[node] for node in nodes], anchor_tips=tips,
        alignment_width=width, anchored_observations=int(states.size),
        unanchored_candidate_residues=len(free_state),
        arrays={key:dict(shape=list(a.shape), dtype=str(a.dtype), sha256=array_digest(a)) for key,a in arrays.items()},
        ancestral_categories_available=True, same_record_sequence_category_correspondence=True,
        legacy_separate_logger_joint_identity_asserted=False, scientific_eligibility=False,
        posterior_qualified=False)
    return summary, arrays


def write_arrays(path, arrays):
    assert set(arrays) == ARRAY_NAMES
    with Path(path).open('xb') as handle:
        np.savez(handle, **arrays)


def verify_arrays(path, expected):
    """Verify all serialized values, dtypes, shapes and keys without pickle."""
    assert set(expected) == ARRAY_NAMES
    with np.load(path, allow_pickle=False) as saved:
        assert len(saved.files) == len(set(saved.files)) and set(saved.files) == ARRAY_NAMES
        for key, value in expected.items():
            actual = saved[key]
            assert actual.dtype == value.dtype and actual.shape == value.shape
            assert np.array_equal(actual, value), key
