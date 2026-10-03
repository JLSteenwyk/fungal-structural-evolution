#!/usr/bin/env python3
"""Independently replay all paired draws, archived native trees and array exports."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from ancestral_chain_attempt import sha
from matched_predictor_branch_fits import config, ROLES
from matched_predictor_branch_inputs import verify
from matched_predictor_resampling import load, add_splits, unpack, SCHEMA, MODES, SUMMARY_FIELDS, read_arrays, values, runtime_caps
from readback_matched_predictor_branch_fits import independent_decode, check_point
from run_matched_predictor_resampling import STATUS as PRODUCER_STATUS

STATUS = 'passed_full_matched_predictor_paired_resampling_independent_readback'


def independent_draw(length, group, mode, replicate):
    digest = hashlib.sha256(f'{SCHEMA}:{group}:{mode}:{replicate}'.encode()).digest()
    rng = np.random.Generator(np.random.PCG64(int.from_bytes(digest[:16], 'big')))
    if mode == 'iid_sites':
        result = rng.integers(length, size=length)
    else:
        assert mode == 'circular_blocks10'
        starts = rng.integers(length, size=(length + 9) // 10)
        result = [(int(start) + offset) % length for start in starts for offset in range(10)][:length]
    return np.asarray(result, dtype='<u2')


def replace_root(value, current, original):
    if isinstance(value, str):
        return original + value[len(current):] if value == current or value.startswith(current + '/') else value
    if isinstance(value, list):
        return [replace_root(item, current, original) for item in value]
    if isinstance(value, dict):
        return {replace_root(key, current, original): replace_root(item, current, original) for key, item in value.items()}
    return value


def case_paths(root, key, mode, replicate):
    folder = root / 'cases' / key / mode
    stem = f'{replicate:04d}'
    return folder / (stem + '.json'), folder / (stem + '.npz'), folder / (stem + '.tar.gz')


def check_manifest(manifest, source):
    expected = [(key, mode, replicate) for key in sorted(source['configs']) for mode in sorted(MODES) for replicate in range(200)]
    assert len(manifest) == len(expected) == 53200
    assert [(row['original_input_id'], row['mode'], row['replicate']) for row in manifest] == expected
    assert len({row['case_id'] for row in manifest}) == len(manifest)
    for row, (key, mode, replicate) in zip(manifest, expected):
        draw = independent_draw(len(source['configs'][key]['columns']), source['resampling_groups'][key], mode, replicate)
        digest = hashlib.sha256(draw.tobytes()).hexdigest()
        identity = hashlib.sha256(json.dumps([SCHEMA, key, mode, replicate, digest], separators=(',', ':')).encode()).hexdigest()
        assert row['case_id'] == identity
        for name, extension in [('receipt', '.json'), ('array', '.npz'), ('archive', '.tar.gz')]:
            assert row[name] == f'cases/{key}/{mode}/{replicate:04d}{extension}'


def read_case(plan, source, key, mode, replicate, root, scratch):
    receipt_path, array_path, archive_path = case_paths(root, key, mode, replicate)
    saved = json.loads(receipt_path.read_text())
    original = source['configs'][key]
    group = source['resampling_groups'][key]
    draw = independent_draw(len(original['columns']), group, mode, replicate)
    draw_hash = hashlib.sha256(draw.tobytes()).hexdigest()
    identity = hashlib.sha256(json.dumps([SCHEMA, key, mode, replicate, draw_hash], separators=(',', ':')).encode()).hexdigest()
    expected_work = str((root / 'work' / identity).resolve())
    expected_seed = int.from_bytes(hashlib.sha256(f'{SCHEMA}:{group}:{mode}:{replicate}'.encode()).digest()[:16], 'big')
    assert saved['schema'] == SCHEMA and saved['case_id'] == identity
    assert (saved['original_input_id'], saved['mode'], saved['replicate']) == (key, mode, replicate)
    assert saved['resampling_group'] == group and saved['seed'] == expected_seed
    assert saved['original_native_work_root'] == expected_work
    assert saved['source_indices_sha256'] == draw_hash
    assert saved['original_config_sha256'] == sha(source['root'] / 'inputs' / key / 'config.json')
    assert saved['array'] == str(array_path) and saved['archive'] == str(archive_path)
    assert sha(array_path) == saved['array_sha256'] and sha(archive_path) == saved['archive_sha256']
    assert saved['scientific_eligibility'] is False
    destination = scratch / identity
    unpack(archive_path, destination, saved['archive_members'])
    assert {p.name for p in destination.iterdir()} == {'inputs', 'native', 'roles'}
    assert {p.name for p in (destination / 'inputs').iterdir()} == {identity}
    assert {p.name for p in (destination / 'native').iterdir()} == {identity}
    folder = destination / 'inputs' / identity
    assert {p.name for p in folder.iterdir()} == {'aa.faa', 'AlphaFold.faa', 'ESMFold.faa', 'topology.nwk', 'config.json'}
    for label, sequences in source['data'][key].items():
        expected = ''.join('>' + taxon + '\n' + ''.join(sequences[taxon][int(i)] for i in draw) + '\n' for taxon in original['taxa'])
        assert (folder / (label + '.faa')).read_text() == expected
    assert (folder / 'topology.nwk').read_bytes() == (source['root'] / 'inputs' / key / 'topology.nwk').read_bytes()
    expected_config = dict(original, input_id=identity,
        columns=[original['columns'][int(i)] for i in draw],
        alignment_sha256={label: sha(folder / (label + '.faa')) for label in source['data'][key]},
        original_input_id=key, resampling_mode=mode, replicate=replicate, resampling_group=group,
        resampling_seed=expected_seed, drawn_source_indices_sha256=draw_hash, schema=SCHEMA)
    expected_config['topology_sha256'] = sha(folder / 'topology.nwk')
    assert json.loads((folder / 'config.json').read_text()) == expected_config
    rows = []
    assert {p.name for p in (destination / 'native' / identity).iterdir()} == {role[0] for role in ROLES}
    assert {p.name for p in (destination / 'roles').iterdir()} == {identity + '-' + role[0] + '.json' for role in ROLES}
    for role in ROLES:
        row = json.loads((destination / 'roles' / (identity + '-' + role[0] + '.json')).read_text())
        native_root = destination / 'native' / identity / role[0]
        assert {p.name for p in native_root.iterdir()} == {'configuration.json', 'attempt.lock', 'attempt-0001'}
        attempt = native_root / 'attempt-0001'
        recipe = config(plan, destination, identity, expected_config, role)
        recipe = replace_root(recipe, str(destination.resolve()), expected_work)
        assert json.loads((native_root / 'configuration.json').read_text()) == recipe
        native = json.loads((attempt / 'receipt.json').read_text())
        assert native['configuration_sha256'] == hashlib.sha256(json.dumps(recipe, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        original_attempt = expected_work + '/native/' + identity + '/' + role[0] + '/attempt-0001'
        command = [value.replace('{attempt}', original_attempt) for value in recipe['command']]
        assert json.loads((attempt / 'command.json').read_text()) == command
        process = json.loads((attempt / 'process.json').read_text())
        assert process['command'] == command
        assert type(process['pid']) is int and process['pid'] > 0 and process['pgid'] == process['pid']
        assert type(process['created']) in (int, float) and process['created'] > 0
        assert {str(p.relative_to(attempt)) for p in attempt.rglob('*') if p.is_file()} == set(native['artifacts']) | {'receipt.json'}
        assert all(sha(attempt / name) == digest for name, digest in native['artifacts'].items())
        assert row['native_receipt'] == original_attempt + '/receipt.json'
        assert row['native_receipt_sha256'] == sha(attempt / 'receipt.json')
        assert row['input_id'] == identity and row['role'] == role[0] and row['seed'] == recipe['seed']
        assert row['source_config_sha256'] == sha(folder / 'config.json')
        assert row['exit_code'] == native['exit_code'] and row['native_status'] == native['status']
        assert row['elapsed_seconds'] == native['elapsed_seconds'] and row['scientific_eligibility'] is False
        if native['exit_code'] != 0:
            assert row['status'] == 'native_unsuccessful_retained' and row['point_estimate'] is None
        elif row['point_estimate'] is not None:
            check_point(row['point_estimate'], independent_decode(attempt, expected_config, source['axes']['positions']))
            assert row['status'] == 'native_point_output_integrity_checked_not_model_qualified'
        else:
            assert row['status'] == 'native_output_requires_review'
        rows.append(row)
    status, lengths = values(rows, source['splits'][key])
    read_arrays(array_path, draw, source['splits'][key], status, lengths)
    counts = dict(Counter(row['status'] for row in rows))
    assert saved['native_role_count'] == 7 and saved['native_status_counts'] == counts
    assert saved['complete'] is bool(np.all(status == 0))
    assert saved['serialized_branch_value_slots'] == lengths.size
    assert saved['finite_branch_values'] == int(np.isfinite(lengths).sum())
    assert saved['native_seconds_sum'] == sum(row['elapsed_seconds'] for row in rows)
    shutil.rmtree(destination)
    return saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    limits = runtime_caps(plan, reader=True)
    source, bindings = load(plan, args.plan)
    add_splits(source, Path(plan['point_output']))
    root = Path(plan['output'])
    producer_path = root / 'receipt.json'
    producer = json.loads(producer_path.read_text())
    assert producer['status'] == PRODUCER_STATUS and producer['plan_sha256'] == sha(args.plan)
    verify(producer['source_hashes'])
    assert not (root / 'readback.json').exists()
    expected_axes = dict(schema=SCHEMA, modes=MODES, replicates_per_mode=200,
        original_inputs=sorted(source['configs']), resampling_groups=source['resampling_groups'],
        native_roles=[role[0] for role in ROLES], splits=source['splits'])
    assert json.loads((root / 'axes.json').read_text()) == expected_axes
    manifest = json.loads((root / 'case_manifest.json').read_text())
    check_manifest(manifest, source)
    expected_files = {root / row[name] for row in manifest for name in ['receipt', 'archive', 'array']}
    assert {p for p in (root / 'cases').rglob('*') if p.is_file()} == expected_files
    scratch = root / 'readback_work'
    scratch.mkdir(exist_ok=False)
    counts = Counter()
    complete = slots = finite = archived = 0
    for index, entry in enumerate(manifest):
        key, mode, replicate = entry['original_input_id'], entry['mode'], entry['replicate']
        paths = case_paths(root, key, mode, replicate)
        for name, path in zip(['receipt', 'array', 'archive'], paths):
            assert str(path.relative_to(root)) == entry[name] and sha(path) == entry[name + '_sha256']
        saved = read_case(plan, source, key, mode, replicate, root, scratch)
        assert entry['case_id'] == saved['case_id']
        counts.update(saved['native_status_counts'])
        complete += saved['complete']; slots += saved['serialized_branch_value_slots']
        finite += saved['finite_branch_values']; archived += len(saved['archive_members'])
        if (index + 1) % 100 == 0:
            print('independent_matched_predictor_resamples_checked', index + 1, '/53200', flush=True)
    assert not list(scratch.iterdir()) and not list((root / 'work').iterdir())
    summary = dict(original_comparison_cases=8750, unique_inputs=133, replicates_per_mode=200, modes=MODES,
        resampling_cases=len(manifest), native_roles=sum(counts.values()), native_status_counts=dict(counts),
        complete_cases=complete, unresolved_cases=len(manifest) - complete,
        serialized_branch_value_slots=slots, finite_branch_values=finite, archived_native_files=archived)
    assert all(producer[key] == summary[key] for key in SUMMARY_FIELDS)
    verify(bindings); verify(producer['source_hashes'])
    assert all(sha(root / name) == digest for name, digest in producer['artifacts'].items())
    bindings[str(producer_path)] = sha(producer_path)
    result = dict(status=STATUS, plan_sha256=sha(args.plan), producer_receipt_sha256=sha(producer_path),
        **summary, source_hashes=bindings, actual_cgroup_limits=limits, scientific_eligibility=False, scope=plan['scope'])
    with (root / 'readback.json').open('x') as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
