"""Decode every saved native sample and verify full residue/count arrays."""
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from independent_native_ancestral_alignment import ALPHABET, fasta_records, tree_labels, native_alignments, project
from run_ortholog_pair_guide_comparison import sha


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def array_digest(value):
    value = np.ascontiguousarray(value, dtype=value.dtype.newbyteorder('<'))
    h = hashlib.sha256(canonical(dict(shape=list(value.shape), dtype=value.dtype.str)))
    h.update(memoryview(value).cast('B')); return h.hexdigest()


def count_states(values, block_coordinates=4096):
    """Fresh full alphabet histogram in bounded coordinate blocks."""
    assert values.ndim == 3 and values.dtype == np.uint8
    draws, nodes, length = values.shape
    assert draws < 65536 and np.all(values < len(ALPHABET))
    coordinates = nodes * length
    temporal = values.reshape(draws, coordinates)
    result = np.empty((coordinates, len(ALPHABET)), dtype=np.uint16)
    for first in range(0, coordinates, block_coordinates):
        last = min(first + block_coordinates, coordinates)
        offsets = np.arange(last-first, dtype=np.int64)[:, None] * len(ALPHABET)
        keys = temporal[:, first:last].T.astype(np.int64) + offsets
        histogram = np.bincount(keys.ravel(), minlength=(last-first) * len(ALPHABET))
        result[first:last] = histogram.reshape(last-first, len(ALPHABET))
    assert np.all(result.sum(axis=1) == draws)
    return result.reshape(nodes, length, len(ALPHABET))


def replay_chain(source, cid, target, readback=False):
    entry = source['details'][cid]
    base = dict(chain_id=cid, model_input_identity=entry['model_input_identity'], seed=entry['seed'], scientific_eligibility=False)
    if entry['status'] == 'unresolved_failed_native_chain_retained':
        return dict(**base, status=entry['status'], native_disposition=entry['native_disposition'], saved_alignments=0, cutoffs={})
    doc = source['doc']
    with Path(entry['input_alignment']).open() as handle: aligned = fasta_records(handle)
    assert len({len(s) for s in aligned.values()}) == 1
    observed = {tip: seq.replace('-', '') for tip, seq in aligned.items()}
    labels, tips = tree_labels(Path(entry['native_files']['runtime-tree.nwk']['path']).read_text())
    assert tips == set(observed) and len(labels) == entry['runtime_nodes'] and len(tips) == entry['tips']
    audit = doc(entry['sample_audit'], entry['sample_audit_sha256'])
    candidates = {}; lengths = {}
    for row in audit['candidate_samples']:
        key = row['iteration'], row['source_node']; assert key not in lengths
        lengths[key] = row['ungapped_length']
        assert row['source_node'] not in candidates.setdefault(row['iteration'], {})
        candidates[row['iteration']][row['source_node']] = row['runtime_node']
    expected_iterations = list(range(0, 1001, 10))
    assert sorted(candidates) == expected_iterations and len(lengths) == 404
    mapping = entry['source_to_runtime_candidates']; assert len(mapping) == len(set(mapping.values())) == 4
    assert all(value == mapping for value in candidates.values())
    assert set(mapping.values()) <= labels - tips
    nodes = sorted(mapping); length = sum(len(seq) for seq in observed.values())
    coords = doc(entry['coordinates']['path'], entry['coordinates']['sha256'])
    assert coords == dict(nodes=nodes, tips=[dict(tip=t,length=len(observed[t])) for t in sorted(observed)],
        alphabet=ALPHABET, input_alignment=entry['input_alignment'],
        input_alignment_sha256=entry['input_alignment_sha256'],
        ordering='state axes: saved iteration, source node, concatenated ungapped input residues in listed tip order; positions are one-based within each tip')
    assert 4 * length == entry['candidate_anchor_coordinates']
    temporary = target.with_suffix(target.suffix + '.partial')
    fresh = np.empty((101, 4, length), dtype=np.uint8)
    unanchored = np.empty((101, 4), dtype=np.int32)
    with np.load(entry['original_state_array']['path'], allow_pickle=False) as stored:
        assert set(stored.files) == {'states', 'iterations', 'unanchored_residue_counts', 'counts_after_250', 'counts_after_500'}
        states, iterations, free = stored['states'], stored['iterations'], stored['unanchored_residue_counts']
        assert states.shape == fresh.shape and states.dtype == np.uint8
        assert iterations.dtype.kind in 'iu' and iterations.tolist() == expected_iterations
        assert free.shape == unanchored.shape and free.dtype == np.int32 and np.all(free >= 0)
        count = 0
        with ExitStack() as stack:
            native = stack.enter_context(Path(entry['native_files']['C1.P1.fastas']['path']).open())
            if readback:
                exported = stack.enter_context(target.open())
            else:
                assert not target.exists()
                target.parent.mkdir(exist_ok=True)
                exported = stack.enter_context(temporary.open('wb'))
            for count, (iteration, sequences) in enumerate(native_alignments(native), 1):
                assert count <= 101 and iteration == expected_iterations[count-1]
                assert set(sequences) == labels, ('Saved runtime node set changed', cid, iteration)
                projected, extra = project(sequences, observed, mapping)
                assert np.array_equal(projected, states[count-1]), ('Projected residues differ', cid, iteration)
                assert np.array_equal(extra, free[count-1]), ('Unanchored counts differ', cid, iteration)
                fresh[count-1] = projected; unanchored[count-1] = extra
                node_lengths = {n:len(sequences[mapping[n]].replace('-', '')) for n in nodes}
                assert all(node_lengths[n] == lengths[(iteration,n)] for n in nodes)
                row = dict(iteration=iteration, alignment_width=len(next(iter(sequences.values()))),
                    runtime_nodes=len(sequences), alignment_sha256=digest(sequences), node_lengths=node_lengths,
                    projected_shape=list(projected.shape), projected_states_sha256=array_digest(projected),
                    unanchored_residue_counts=extra.tolist(), scientific_eligibility=False)
                if readback:
                    line = exported.readline(); assert line, ('Missing exported native sample', cid, iteration)
                    assert json.loads(line) == row, ('Changed exported native sample', cid, iteration)
                else: exported.write(canonical(row) + b'\n')
            assert count == 101
            if readback: assert exported.readline() == '', ('Extra exported native sample', cid)
            else: exported.flush(); os.fsync(exported.fileno())
        cutoffs = {}
        for cut, draws in [('250',75),('500',50)]:
            counts = count_states(fresh[np.asarray(expected_iterations) > int(cut)])
            original = stored['counts_after_' + cut]
            assert original.dtype == np.uint16 and original.shape == counts.shape
            assert np.array_equal(original, counts), ('Cutoff count arrays differ', cid, cut)
            cutoffs[cut] = dict(retained_samples=draws, count_shape=list(counts.shape), count_cells=int(counts.size),
                counts_sha256=array_digest(counts), state_totals=counts.sum(axis=(0,1),dtype=np.uint64).tolist())
    assert int(unanchored.sum()) == entry['inherited_unanchored_residue_observations']
    if not readback: os.replace(temporary, target)
    return dict(**base, status='all_native_saved_alignments_and_projections_replayed_not_posterior_qualification',
        sample_audit_sha256=entry['sample_audit_sha256'], input_alignment_sha256=entry['input_alignment_sha256'],
        native_files=entry['native_files'], source_state_array=entry['original_state_array'],
        source_to_runtime_candidates=mapping, saved_alignments=101, candidate_frames=404,
        state_observations=int(fresh.size), candidate_anchor_coordinates=4*length,
        unanchored_residue_observations=int(unanchored.sum()), cutoff_count_cells=sum(c['count_cells'] for c in cutoffs.values()),
        full_projected_states_sha256=array_digest(fresh), full_unanchored_counts_sha256=array_digest(unanchored),
        frames_path=str(target.name), frames_sha256=sha(target), cutoffs=cutoffs,
        independent_source_to_runtime_clade_mapping=False)
