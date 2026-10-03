"""Manual replay of every short-sampler alignment and available site properties.

Ancestral category labels are absent from the native logger. No complete
ancestral category trajectory, posterior adequacy or likelihood proof is claimed.
"""
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from independent_native_ancestral_alignment import fasta_records, native_alignments, project
from independent_native_ancestral_topology import match_trees


NATIVE_ALPHABET = 'ARNDCQEGHILKMFPSTWYV'
SAVED_ITERATIONS = [0, 10, 20]
SUCCESS = 'all_short_native_alignments_and_available_site_properties_replayed_not_posterior'
SUMMARY_FIELDS = ['full_chains', 'full_groups', 'effective_inputs', 'original_configuration_aliases',
    'checked_chains', 'unresolved_chains', 'complete_groups', 'unresolved_groups',
    'intact_chains_in_unresolved_groups', 'saved_alignments', 'candidate_frames',
    'projected_state_observations', 'tip_category_state_pairs', 'rate_property_cells',
    'unanchored_candidate_residue_observations', 'node_pairs', 'nonroot_branch_pairs',
    'negative_branch_pairs_require_review', 'status_counts', 'ancestral_categories_available',
    'posterior_qualified']


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def array_digest(value):
    return hashlib.sha256(value.tobytes(order='C')).hexdigest()


def strict_json(text):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result: raise ValueError('Duplicate JSON key: ' + key)
            result[key] = value
        return result
    return json.loads(text, object_pairs_hook=pairs,
                      parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)))


def site_frame(frame, iteration, sequences, observed):
    assert set(frame) == {'iter', 'catStates', 'properties', 'conditions'}
    assert type(frame['iter']) is int and frame['iter'] == iteration
    assert set(frame['catStates']) == set(observed), 'Missing, extra or ancestral category records'
    assert frame['conditions'] == {} and set(frame['properties']) == {'rate'}
    rates = frame['properties']['rate']
    assert len(rates) == 4 and all(len(row) == 20 for row in rates)
    for row in rates:
        assert all(type(x) in [int, float] and math.isfinite(x) and x >= 0 for x in row)
        assert all(x == row[0] for x in row), 'Rate unexpectedly depends on amino-acid state'
    assert math.isclose(sum(row[0] for row in rates) / 4, 1, rel_tol=1e-10, abs_tol=1e-10)
    pairs = 0
    for tip in sorted(observed):
        value = frame['catStates'][tip]
        assert set(value) == {'categories', 'states'}
        letters = sequences[tip].replace('-', '')
        categories, states = value['categories'], value['states']
        assert len(categories) == len(states) == len(letters) == len(observed[tip])
        assert all(type(x) is int and 0 <= x < 4 for x in categories)
        assert all(type(x) is int and 0 <= x < 20 for x in states)
        assert ''.join(NATIVE_ALPHABET[x] for x in states) == letters, 'Tip state/alignment disagreement'
        pairs += len(states)
    return dict(tip_category_state_pairs=pairs, rate_property_cells=80,
                complete_site_property_record_sha256=digest(frame),
                category_rates=[row[0] for row in rates], ancestral_categories_available=False)


def mapping_rows(path, chain):
    dataset = 'whole' if '-whole-' in chain['original_configuration_ids'][0] else 'domain'
    with Path(path).open() as handle:
        rows = [r for r in csv.DictReader(handle, delimiter='\t') if
                (r['guide'], r['family'], r['dataset']) == ('profile', chain['family'], dataset)]
    assert len(rows) == 4 and {int(r['level']) for r in rows} == set(range(4))
    return sorted(rows, key=lambda r: int(r['level']))


def replay(job, disposition, mapping, array_path=None):
    """Return independently decoded records; optionally write fresh projection arrays."""
    chain = job['chain']; cid = chain['chain_id']
    expected = dict(chain_id=cid, effective_input_group=chain['effective_input_group'],
        model_input_identity=chain['effective_input_group'] + '-' + chain['prior_label'],
        prior_label=chain['prior_label'], chain_role=chain['chain'], family=chain['family'],
        original_configuration_ids=chain['original_configuration_ids'], seed=chain['seed'],
        source_seed=job['source_seed'], scientific_eligibility=False, posterior_qualified=False)
    assert all(disposition[k] == value for k, value in expected.items())
    receipt_path = Path(disposition['native_receipt'])
    assert sha(receipt_path) == disposition['native_receipt_sha256']
    receipt = strict_json(receipt_path.read_text()); base = receipt_path.parent
    assert strict_json((base.parent / 'configuration.json').read_text()) == job['config']
    assert strict_json((base / 'command.json').read_text()) == job['config']['command']
    assert strict_json((base / 'process.json').read_text())['command'] == job['config']['command']
    assert receipt['configuration_sha256'] == digest(job['config'])
    for name, h in receipt['artifacts'].items(): assert sha(base / name) == h
    result = dict(**expected, native_receipt=str(receipt_path), native_receipt_sha256=sha(receipt_path),
        inherited_native_disposition=disposition['status'], native_exit_code=receipt['exit_code'],
        ancestral_categories_available=False)
    if disposition['status'] != 'full_short_sampler_output_integrity_checked_not_posterior':
        assert disposition['status'] in ['unsuccessful_sampler_qualification_attempt_retained',
                                         'invalid_sampler_qualification_output_retained']
        return dict(**result, status='unresolved_original_short_sampler_disposition_retained',
                    frames=[], node_rows=[], candidate_mapping={}, projected_state_observations=0,
                    unanchored_candidate_residue_observations=0, negative_branch_pairs_require_review=0)
    assert receipt['exit_code'] == 0
    directories = list(base.glob('independent-chain-*')); assert len(directories) == 1
    directory = directories[0]
    for name in ['tree', 'alignment']: assert sha(chain[name]) == chain[name + '_sha256']
    matched = match_trees(Path(chain['tree']).read_text(), (directory / 'runtime-tree.nwk').read_text())
    with Path(chain['alignment']).open() as handle:
        observed = {name: sequence.replace('-', '') for name, sequence in fasta_records(handle).items()}
    assert sorted(observed) == matched['tips'] and len(observed) == chain['proteins']
    bits = {tip: 1 << i for i, tip in enumerate(matched['tips'])}; candidates = {}; levels = {}
    for row in mapping_rows(mapping, chain):
        retained = strict_json(row['retained_set_json'])
        assert retained and len(retained) == len(set(retained)) and set(retained) <= set(bits)
        if 'retained_descendants' in row: assert len(retained) == int(row['retained_descendants'])
        mask = sum(bits[tip] for tip in retained); node = row['source_node']; level = row['level']
        assert node not in candidates and matched['source_labels'][node]['mask'] == mask
        assert (int(level) == 3) == (mask == (1 << len(bits)) - 1)
        candidates[node] = matched['runtime_index'][mask]['label']; levels[node] = level
    assert len(candidates) == len(set(candidates.values())) == 4
    audit = disposition['sample_audit']; assert audit['mapping_sha256'] == sha(mapping)
    lengths = {}
    for row in audit['candidate_samples']:
        key = row['iteration'], row['source_node']; assert key not in lengths
        assert row['runtime_node'] == candidates[row['source_node']] and row['level'] == levels[row['source_node']]
        lengths[key] = row['ungapped_length']
    assert set(lengths) == {(i, node) for i in SAVED_ITERATIONS for node in candidates}
    with (directory / 'C1.log').open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        assert len(reader.fieldnames) == len(set(reader.fieldnames))
        logs = list(reader)
    assert [int(row['iter']) for row in logs] == list(range(21))
    assert all(None not in row and all(value is not None for value in row.values()) for row in logs)
    for row in logs:
        scores = [float(row[name]) for name in ['prior', 'likelihood', 'posterior']]
        assert all(math.isfinite(v) for v in scores) and abs(scores[0] + scores[1] - scores[2]) < 1e-7
    frames = []; states = []; free = []
    with (directory / 'C1.P1.fastas').open() as alignments, (directory / 'C1.P1.site-property-samples.jsonl').open() as properties:
        for index, (iteration, sequences) in enumerate(native_alignments(alignments)):
            assert index < 3 and iteration == SAVED_ITERATIONS[index]
            assert set(sequences) == set(matched['runtime_labels'])
            projected, unanchored = project(sequences, observed, candidates)
            node_lengths = {node: len(sequences[target].replace('-', '')) for node, target in candidates.items()}
            assert all(node_lengths[node] == lengths[(iteration, node)] for node in candidates)
            line = properties.readline(); assert line, 'Missing site-property sample'
            property_record = site_frame(strict_json(line), iteration, sequences, observed)
            frames.append(dict(iteration=iteration, alignment_width=len(next(iter(sequences.values()))),
                alignment_sha256=digest(sequences), node_lengths=node_lengths,
                projected_shape=list(projected.shape), projected_states_sha256=array_digest(projected),
                unanchored_residue_counts=unanchored.tolist(), **property_record, scientific_eligibility=False))
            states.append(projected); free.append(unanchored)
        assert len(frames) == 3 and properties.readline() == '', 'Missing or extra saved sample'
    state_array, free_array = np.stack(states), np.stack(free)
    if array_path is not None:
        with Path(array_path).open('xb') as handle:
            np.savez(handle, states=state_array, iterations=np.array(SAVED_ITERATIONS, dtype=np.int32),
                     unanchored_residue_counts=free_array)
    result.update(status=SUCCESS, frames=frames, node_rows=matched['rows'], candidate_mapping=candidates,
        projected_state_observations=int(state_array.size), full_projected_states_sha256=array_digest(state_array),
        full_unanchored_counts_sha256=array_digest(free_array),
        unanchored_candidate_residue_observations=int(free_array.sum()),
        negative_branch_pairs_require_review=matched['negative_branch_pairs_require_review'])
    return result, state_array, free_array


def summarize(rows, jobs):
    expected = {j['chain']['chain_id']: j for j in jobs}
    assert len(rows) == len(expected) == 1620 and len({r['chain_id'] for r in rows}) == 1620
    assert set(expected) == {r['chain_id'] for r in rows}
    groups = defaultdict(list); counts = Counter(); inputs = set(); aliases = set()
    for row in rows:
        chain = expected[row['chain_id']]['chain']
        assert row['seed'] == chain['seed'] and row['scientific_eligibility'] is False
        assert row['posterior_qualified'] is row['ancestral_categories_available'] is False
        assert row['model_input_identity'] == chain['effective_input_group'] + '-' + chain['prior_label']
        assert row['original_configuration_ids'] == chain['original_configuration_ids']
        groups[row['model_input_identity']].append(row); inputs.add(chain['effective_input_group'])
        aliases.update(chain['original_configuration_ids']); counts['status:' + row['status']] += 1
        if row['status'] != SUCCESS:
            assert row['status'] == 'unresolved_original_short_sampler_disposition_retained'
            assert row['frames'] == row['node_rows'] == [] and row['candidate_mapping'] == {}
            counts['unresolved_chains'] += 1; continue
        assert [r['iteration'] for r in row['frames']] == SAVED_ITERATIONS
        counts['checked_chains'] += 1; counts['saved_alignments'] += 3; counts['candidate_frames'] += 12
        for key in ['projected_state_observations', 'unanchored_candidate_residue_observations',
                    'negative_branch_pairs_require_review']: counts[key] += row[key]
        counts['node_pairs'] += len(row['node_rows'])
        counts['nonroot_branch_pairs'] += sum(not node['is_root'] for node in row['node_rows'])
        for frame in row['frames']:
            for key in ['tip_category_state_pairs', 'rate_property_cells']: counts[key] += frame[key]
    assert len(groups) == 405 and len(inputs) == 135 and len(aliases) == 324
    for group in groups.values():
        assert len(group) == 4 and {r['chain_role'] for r in group} == {1, 2, 3, 4}
        complete = all(r['status'] == SUCCESS for r in group)
        counts['complete_groups' if complete else 'unresolved_groups'] += 1
        if not complete: counts['intact_chains_in_unresolved_groups'] += sum(r['status'] == SUCCESS for r in group)
    result = {key: int(counts[key]) for key in SUMMARY_FIELDS if key not in
              ['full_chains', 'full_groups', 'effective_inputs', 'original_configuration_aliases',
               'status_counts', 'ancestral_categories_available', 'posterior_qualified']}
    return dict(result, full_chains=1620, full_groups=405, effective_inputs=135, original_configuration_aliases=324,
        status_counts=dict(Counter(row['status'] for row in rows)), ancestral_categories_available=False,
        posterior_qualified=False)
