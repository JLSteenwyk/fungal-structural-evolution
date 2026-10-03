#!/usr/bin/env python3
"""Qualify every full-grid draw and native/archive/serialization software contracts."""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import tarfile

import numpy as np

from ancestral_chain_attempt import sha
from matched_predictor_branch_inputs import verify
from matched_predictor_resampling import load, add_splits, indices, case_id, perform, values, read_arrays, unpack, SCHEMA, MODES
from readback_matched_predictor_branch_inputs import raw_fasta
from readback_matched_predictor_resampling import independent_draw, read_case, check_manifest


def reject(action):
    try:
        action()
    except (AssertionError, ValueError, KeyError, FileNotFoundError, OSError):
        return
    raise AssertionError('Malformed resampling contract accepted')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists() and not args.receipt.exists()
    args.output.mkdir(parents=True)
    point_plan_path = Path('metadata/matched_predictor_branch_fits_plan_20261003_v1.json')
    point_plan = json.loads(point_plan_path.read_text())
    own = ['matched_predictor_resampling', 'run_matched_predictor_resampling',
           'readback_matched_predictor_resampling', 'check_matched_predictor_resampling']
    pins = {f'scripts/{name}.py': sha(f'scripts/{name}.py') for name in own}
    source_plan = dict(point_plan, pins=pins, point_completion='metadata/matched_predictor_branch_fits_completed_20261003_v1.json',
        point_output=point_plan['output'], replicates_per_mode=200, modes=MODES)
    path = args.output / 'source_plan.json'
    path.write_text(json.dumps(source_plan, indent=2) + '\n')
    source, bindings = load(source_plan, path)
    add_splits(source, Path(point_plan['output']))
    manifest = []
    total_slots = 0
    for key in sorted(source['configs']):
        cfg = source['configs'][key]; group = source['resampling_groups'][key]
        for mode in sorted(MODES):
            for replicate in range(200):
                draw = indices(len(cfg['columns']), group, mode, replicate)
                assert np.array_equal(draw, independent_draw(len(cfg['columns']), group, mode, replicate))
                assert draw.dtype == np.dtype('<u2') and draw.shape == (len(cfg['columns']),)
                assert int(draw.max()) < len(cfg['columns'])
                identity = case_id(key, mode, replicate, draw)
                base = f'cases/{key}/{mode}/{replicate:04d}'
                manifest.append(dict(original_input_id=key, mode=mode, replicate=replicate, case_id=identity,
                    receipt=base + '.json', array=base + '.npz', archive=base + '.tar.gz'))
                total_slots += 7 * len(source['splits'][key])
    check_manifest(manifest, source)
    serialized = args.output / 'full_draw_grid.json'
    serialized.write_text(json.dumps(manifest, indent=2) + '\n')
    check_manifest(json.loads(serialized.read_text()), source)
    negatives = []
    for label in ['missing_case', 'duplicate_case', 'changed_identity', 'changed_mode', 'foreign_archive_path']:
        bad = copy.deepcopy(manifest)
        if label == 'missing_case': bad.pop()
        elif label == 'duplicate_case': bad[-1] = bad[0]
        elif label == 'changed_identity': bad[0]['case_id'] = '0' * 64
        elif label == 'changed_mode': bad[0]['mode'] = 'unpaired'
        else: bad[0]['archive'] = '../foreign.tar.gz'
        reject(lambda: check_manifest(bad, source)); negatives.append(label)
    prior_gate_path = Path('metadata/matched_predictor_branch_fits_software_validation_20261003_v2.json')
    prior_gate = json.loads(prior_gate_path.read_text())
    fixture_config_path = Path(next(name for name in prior_gate['artifacts'] if '/fixture_inputs/' in name and name.endswith('/config.json')))
    assert sha(fixture_config_path) == prior_gate['artifacts'][str(fixture_config_path)]
    fixture_config = json.loads(fixture_config_path.read_text())
    key = fixture_config['input_id']; fixture_root = fixture_config_path.parent.parent.parent
    fixture = dict(root=fixture_root, configs={key: fixture_config},
        data={key: {label: raw_fasta(fixture_config_path.parent / (label + '.faa')) for label in ['aa', 'AlphaFold', 'ESMFold']}},
        resampling_groups={key: hashlib.sha256(b'paired-resampling-synthetic-software-v1').hexdigest()},
        splits={key: ['0x1', '0x2', '0x3', '0x4', '0x8', '0x10', '0x18']},
        axes={'positions': {taxon: i for i, taxon in enumerate(fixture_config['taxa'])}})
    fixture_plan = dict(executable=point_plan['executable'], models=point_plan['models'], resources=point_plan['resources'])
    native_root = args.output / 'native_fixtures'; (native_root / 'work').mkdir(parents=True)
    scratch = args.output / 'readback_work'; scratch.mkdir()
    cases = []
    for mode in MODES:
        saved = perform(fixture_plan, fixture, key, mode, 0, native_root)
        assert saved['complete'] is True
        assert read_case(fixture_plan, fixture, key, mode, 0, native_root, scratch) == saved
        cases.append(saved)
    forced_plan = copy.deepcopy(fixture_plan)
    forced_plan['resources']['native_cpu_seconds'] = 0
    saved = perform(forced_plan, fixture, key, 'iid_sites', 1, native_root)
    assert saved['complete'] is False and saved['finite_branch_values'] == 0
    assert saved['native_status_counts'] == {'native_unsuccessful_retained': 7}
    assert read_case(forced_plan, fixture, key, 'iid_sites', 1, native_root, scratch) == saved
    cases.append(saved)
    # Explicit review and failed dispositions keep NaN slots instead of dropping rows.
    archive_rows = []
    with tarfile.open(cases[0]['archive'], 'r:gz') as archive:
        for member in archive:
            if member.name.startswith('roles/'):
                archive_rows.append(json.loads(archive.extractfile(member).read()))
    mixed = copy.deepcopy(archive_rows)
    mixed[0].update(status='native_unsuccessful_retained', point_estimate=None)
    mixed[1].update(status='native_output_requires_review', point_estimate=None)
    status, lengths = values(mixed, fixture['splits'][key])
    assert sorted(status.tolist()) == [0, 0, 0, 0, 0, 1, 2]
    draw = indices(240, fixture['resampling_groups'][key], MODES[0], 0)
    mixed_path = args.output / 'mixed_failed_review_arrays.npz'
    with mixed_path.open('xb') as handle:
        np.savez(handle, source_indices=draw, role_status=status, branch_lengths=lengths)
    read_arrays(mixed_path, draw, fixture['splits'][key], status, lengths)
    for label in ['changed_indices', 'wrong_dtype', 'missing_array', 'negative_branch', 'failed_finite_branch', 'changed_status']:
        bad = dict(source_indices=draw.copy(), role_status=status.copy(), branch_lengths=lengths.copy())
        if label == 'changed_indices': bad['source_indices'][0] = (int(draw[0]) + 1) % 240
        elif label == 'wrong_dtype': bad['source_indices'] = bad['source_indices'].astype('int64')
        elif label == 'missing_array': del bad['role_status']
        elif label == 'negative_branch': bad['branch_lengths'][np.flatnonzero(status == 0)[0], 0] = -1
        elif label == 'failed_finite_branch': bad['branch_lengths'][np.flatnonzero(status != 0)[0], 0] = 0
        else: bad['role_status'][0] = 3
        p = args.output / (label + '.npz')
        with p.open('xb') as handle: np.savez(handle, **bad)
        reject(lambda: read_arrays(p, draw, fixture['splits'][key], status, lengths)); negatives.append(label)
    for label in ['parent_path', 'absolute_path', 'symlink', 'wrong_member_hash', 'missing_member']:
        p = args.output / (label + '.tar.gz')
        name = '../escape' if label == 'parent_path' else '/escape' if label == 'absolute_path' else 'file'
        data = b'checked software fixture'
        with tarfile.open(p, 'w:gz') as archive:
            item = tarfile.TarInfo(name); item.size = len(data)
            if label == 'symlink': item.type = tarfile.SYMTYPE; item.linkname = '../escape'
            archive.addfile(item, io.BytesIO(data))
        expected = {name: hashlib.sha256(data).hexdigest()}
        if label == 'wrong_member_hash': expected[name] = '0' * 64
        if label == 'missing_member': expected['extra'] = hashlib.sha256(data).hexdigest()
        reject(lambda: unpack(p, args.output / (label + '_extract'), expected)); negatives.append(label)
    assert not list((native_root / 'work').iterdir()) and not list(scratch.iterdir())
    bindings[str(prior_gate_path)] = sha(prior_gate_path)
    bindings.update(prior_gate['artifacts']); verify(bindings)
    result = dict(status='passed_full_matched_predictor_paired_resampling_software_contracts',
        checked_utc=datetime.now(timezone.utc).isoformat(), full_real_inputs=133, original_comparison_cases=8750,
        full_real_resampling_cases=len(manifest), future_native_roles=7 * len(manifest),
        serialized_branch_value_slots=total_slots, distinct_alignment_draw_groups=len(set(source['resampling_groups'].values())),
        modes=MODES, replicates_per_mode=200, full_draw_serialization_checked=True,
        synthetic_native_cases=3, synthetic_native_roles=21, synthetic_unsuccessful_roles_retained=7,
        artificial_failed_and_review_array_dispositions=2, malformed_cases_rejected=negatives,
        source_hashes=bindings, artifacts={str(p): sha(p) for p in args.output.rglob('*') if p.is_file()},
        scientific_eligibility=False,
        scope='All 53,200 real draw identities and independent site/block indices checked and serialized, retaining the closed full 8,750-case input scope. Two paired synthetic native cases (fourteen fits) and one software-only CPU-zero failure case (seven fits) exercise lossless native archives, independent DendroPy tree/report readback and exact array exports. Sixteen malformed contracts rejected. Fixed predictions/topologies and selected overlap remain conditioning inputs. No biological pilot, accepted evolutionary effects, posterior qualification or GPU work.')
    with args.receipt.open('x') as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ['source_hashes', 'artifacts']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
