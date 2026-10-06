#!/usr/bin/env python3
"""Replay the full derived census using NumPy values and separate integer costs."""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind


KEYS = {'states', 'categories', 'anchor_columns', 'anchor_tip_indices', 'anchor_tip_offsets',
        'unanchored_candidate_indices', 'unanchored_columns', 'unanchored_states',
        'unanchored_categories', 'category_rates'}
SUCCESS = 'full_joint_short_sampler_output_integrity_checked_not_posterior'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('producer', 'producer-transport', 'receipt'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    producer = json.loads(args.producer.read_text())
    transport = json.loads(args.producer_transport.read_text())
    assert producer['status'] == 'complete_full_ancestral_derived_storage_census_pending_readback'
    assert transport['original_tool_terminal_exit_code'] == 0
    assert transport['validation_sha256'] == sha(args.producer)
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    pins = dict(transport['source_hashes'])
    tables = {}
    for name in ('attempt_table', 'file_table', 'member_table', 'derived_file_table'):
        with Path(producer[name]).open() as handle:
            tables[name] = list(csv.DictReader(handle, delimiter='\t'))
    attempt_table = {(r['phase'], r['chain_id']): r for r in tables['attempt_table']}
    file_table = {r['path']: r for r in tables['file_table']}
    member_table = {(r['path'], r['member']): r for r in tables['member_table']}
    leaf_table = {r['path']: r for r in tables['derived_file_table']}
    assert len(attempt_table) == len(tables['attempt_table']) == 1644
    assert len(file_table) == len(tables['file_table']) == 4818
    assert len(member_table) == len(tables['member_table']) == 48180
    assert len(leaf_table) == len(tables['derived_file_table']) == 6477
    roots = []
    source = {}
    failed = set()
    jobs10 = {}
    for phase, name in [('V6_original', 'metadata/baliphy_scalar_v6_sampler_execution_plan_20261004_v1.json'),
                        ('V10_comparison', 'metadata/baliphy_log_alpha_v10_full_grid_plan_20261005_v1.json')]:
        assert sha(name) == pins[name]
        plan = json.loads(Path(name).read_text())
        root = Path(plan['output'])
        roots.append((phase, root))
        rows = json.loads((root / 'dispositions.json').read_text())
        assert sha(root / 'dispositions.json') == pins[str(root / 'dispositions.json')]
        for row in rows:
            source[phase, row['chain_id']] = row
            assert json.loads((root / 'chains' / (row['chain_id'] + '.json')).read_text()) == row
        if phase == 'V6_original':
            failed = {r['chain_id'] for r in rows if r['exit_code'] != 0}
        else:
            assert sha(plan['jobs']) == pins[plan['jobs']]
            jobs10 = {j['chain']['chain_id']: j for j in json.loads(Path(plan['jobs']).read_text())}
    assert len(failed) == 24 and {j['source_v6_chain_id'] for j in jobs10.values()} == failed
    assert set(source) == set(attempt_table)
    leaf_paths, expected_array_paths = set(), set()
    totals = Counter()
    measured_role_files = defaultdict(list)
    raw_bytes = overhead_bytes = 0
    for phase, root in roots:
        assert {p.name for p in root.iterdir() if p.is_dir()} == {'attempts', 'chains', 'frames'}
        leaves = [p for p in root.iterdir() if p.is_file()]
        leaves += list((root / 'chains').glob('*.json'))
        leaves += list((root / 'frames').glob('*/*.npz'))
        assert {p.name for p in (root / 'chains').iterdir()} == {cid + '.json' for ph, cid in source if ph == phase}
        assert {p.name for p in (root / 'frames').iterdir()} == {cid for (ph, cid), r in source.items() if ph == phase and r['status'] == SUCCESS}
        for p in leaves:
            assert str(p) not in leaf_paths
            leaf_paths.add(str(p))
            r = leaf_table[str(p)]
            assert r['phase'] == phase and int(r['bytes']) == p.stat().st_size
            assert r['sha256'] == sha(p) == pins[str(p)]
            kind = 'projection_array' if p.suffix == '.npz' else 'role_disposition' if p.parent.name == 'chains' else 'stage_metadata'
            assert r['kind'] == kind
            if p.suffix == '.lock':
                assert p.stat().st_size == 0 and r['custody'] == 'observed_zero_byte_lock_not_completion_evidence'
            elif phase == 'V10_comparison' and p == root / 'receipt.json':
                assert r['custody'] == 'newly_measured_supplementary_stage_metadata_not_completion_evidence'
            elif p.name == 'completion_archive.json':
                assert r['custody'] in {'closed_completion_archive', 'closed_original_source'}
                completion = json.loads(Path('metadata/baliphy_scalar_v6_short_sampler_completed_20261004_v1.json').read_text())
                assert r['sha256'] == completion['full_hash_archive_sha256']
            else:
                assert r['custody'] == 'closed_original_source'
            totals[kind] += p.stat().st_size
    assert leaf_paths == set(leaf_table)
    for (phase, cid), row in source.items():
        ref = attempt_table[phase, cid]
        baseline = phase == 'V10_comparison' or cid not in failed
        original = cid if phase == 'V6_original' else jobs10[cid]['source_v6_chain_id']
        for name in ('family', 'effective_input_group', 'prior_label', 'chain_role'):
            assert ref[name] == str(row[name])
        assert ref['original_v6_chain_id'] == original
        assert ref['baseline_for_cost_scenario'] == str(baseline)
        assert ref['original_integrity_disposition'] == row['status'] and int(ref['native_exit_code']) == row['exit_code']
        assert ref['scientifically_imputed'] == ref['scientific_eligibility'] == ref['posterior_qualified'] == 'False'
        assert row['scientific_eligibility'] is row['posterior_qualified'] is False
        frames = row['joint_frames']
        if row['status'] == SUCCESS:
            assert [f['iteration'] for f in frames] == [0, 10, 20] and row['exit_code'] == 0
            folder = Path(frames[0]['projection_array']).parent
            assert {p.name for p in folder.iterdir()} == {'frame-0.npz', 'frame-10.npz', 'frame-20.npz'}
        else:
            assert not frames
        for frame in frames:
            p = frame['projection_array']
            assert p not in expected_array_paths
            expected_array_paths.add(p)
            entry = file_table[p]
            assert entry['phase'] == phase and entry['chain_id'] == cid and int(entry['iteration']) == frame['iteration']
            assert entry['sha256'] == frame['projection_array_sha256'] == pins[p]
            assert int(entry['bytes']) == Path(p).stat().st_size
            assert entry['scientific_eligibility'] == entry['posterior_qualified'] == 'False'
            data_size, headers = 0, 0
            with np.load(p, allow_pickle=False) as saved:
                assert len(saved.files) == 10 and set(saved.files) == KEYS == set(frame['arrays'])
                for key in sorted(KEYS):
                    array = saved[key]
                    wanted = frame['arrays'][key]
                    member = member_table[p, key]
                    assert array.flags.c_contiguous and not array.dtype.hasobject
                    assert list(array.shape) == wanted['shape'] == json.loads(member['shape'])
                    assert str(array.dtype) == wanted['dtype'] == member['dtype']
                    digest = hashlib.sha256(memoryview(array).cast('B')).hexdigest() if array.size else hashlib.sha256(b'').hexdigest()
                    assert digest == wanted['sha256'] == member['array_data_sha256']
                    assert array.nbytes == math.prod(array.shape) * array.dtype.itemsize == int(member['raw_array_bytes'])
                    assert member['phase'] == phase and member['chain_id'] == cid and int(member['iteration']) == frame['iteration']
                    assert member['scientific_eligibility'] == member['posterior_qualified'] == 'False'
                    # Independent NPY header measurement: read the public version/length prefix directly.
                    with saved.zip.open(key + '.npy') as handle:
                        prefix = handle.read(10)
                        assert prefix[:8] == b'\x93NUMPY\x01\x00'
                        header = 10 + int.from_bytes(prefix[8:10], 'little')
                        info = saved.zip.getinfo(key + '.npy')
                        assert info.compress_type == 0 and info.file_size == array.nbytes + header
                    assert header == int(member['npy_header_bytes']) and info.file_size == int(member['member_bytes'])
                    data_size += array.nbytes
                    headers += header
            zip_bytes = int(entry['bytes']) - data_size - headers
            assert data_size == int(entry['raw_array_bytes']) and headers == int(entry['npy_header_bytes'])
            assert zip_bytes == int(entry['zip_overhead_bytes']) and zip_bytes >= 0
            raw_bytes += data_size
            overhead_bytes += headers + zip_bytes
            measured_role_files[phase, cid].append(entry)
        measured = measured_role_files[phase, cid]
        for name, value in [('measured_frames', len(measured)),
                            ('measured_array_file_bytes', sum(int(r['bytes']) for r in measured)),
                            ('measured_raw_array_bytes', sum(int(r['raw_array_bytes']) for r in measured)),
                            ('measured_serialization_overhead_bytes', sum(int(r['npy_header_bytes']) + int(r['zip_overhead_bytes']) for r in measured)),
                            ('maximum_observed_frame_bytes', max((int(r['bytes']) for r in measured), default=0))]:
            assert int(ref[name]) == value
    assert expected_array_paths == set(file_table)
    assert set(member_table) == {(p, k) for p in file_table for k in KEYS}
    baselines = [r for r in attempt_table.values() if r['baseline_for_cost_scenario'] == 'True']
    assert len(baselines) == len({r['original_v6_chain_id'] for r in baselines}) == 1620
    assert len({r['effective_input_group'] for r in baselines}) == 135
    assert len({(r['effective_input_group'], r['prior_label']) for r in baselines}) == 405
    assert Counter(r['prior_label'] for r in baselines) == {'broad': 540, 'centered': 540, 'package': 540}
    own = proxy = 0
    proxy_rows = []
    for key, row in attempt_table.items():
        if row['baseline_for_cost_scenario'] == 'False':
            assert int(row['projected_array_bytes']) == 0 and row['projection_basis'] == 'original_native_failure_retained_not_baseline'
        elif int(row['measured_frames']):
            quotient, remainder = divmod(int(row['measured_array_file_bytes']) * 1001, 3)
            value = quotient + (1 if remainder else 0)
            assert value == int(row['projected_array_bytes'])
            assert row['projection_basis'] == 'own_three_integrity_checked_frames'
            assert row['resource_proxy_donor_chain_id'] == '' and int(row['resource_proxy_frame_bytes']) == 0
            own += value
        else:
            candidates = [r for r in baselines if r['effective_input_group'] == row['effective_input_group']
                          and r['prior_label'] == row['prior_label'] and int(r['measured_frames'])]
            ranked = sorted(candidates, key=lambda r: (int(r['maximum_observed_frame_bytes']), r['chain_id']))
            donor = ranked[-1]
            assert row['projection_basis'] == 'same_input_prior_peer_maximum_frame_resource_proxy_only'
            assert row['resource_proxy_donor_chain_id'] == donor['chain_id']
            size = int(donor['maximum_observed_frame_bytes'])
            assert int(row['resource_proxy_frame_bytes']) == size
            assert int(row['projected_array_bytes']) == size * 1001
            proxy += size * 1001
            proxy_rows.append(row['chain_id'])
    assert len(proxy_rows) == 14 and sum(int(r['measured_frames']) > 0 for r in baselines) == 1606
    native_path = 'metadata/current_ancestral_horizon_cost_20261006_v1.json'
    assert sha(native_path) == pins[native_path]
    native = json.loads(Path(native_path).read_text())['linear_native_total_with_latent_allowance_bytes']
    checks = dict(measured_attempts=1644, original_roles=1620, effective_inputs=135, original_quartets=405,
        failed_originals_retained=24, comparison_attempts_retained=24, computational_baselines=1620,
        integrity_checked_baselines=1606, review_baselines_without_exports=14, measured_projection_files=4818,
        measured_array_members=48180, measured_derived_files=6477, measured_bytes_by_kind=dict(totals),
        measured_derived_bytes=sum(totals.values()), measured_raw_array_bytes=raw_bytes,
        measured_serialization_overhead_bytes=overhead_bytes,
        largest_observed_projection_file_bytes=max(int(r['bytes']) for r in file_table.values()),
        proposed_iterations=10000, saved_frame_interval=10, proposed_frames_per_role=1001,
        linear_observed_role_array_projection_bytes=own, resource_only_peer_proxy_projection_bytes=proxy,
        array_projection_with_review_resource_proxies_bytes=own + proxy,
        native_projection_with_latent_allowance_bytes=native,
        native_plus_array_projection_with_resource_proxies_bytes=native + own + proxy)
    for k, value in checks.items():
        assert producer[k] == value, k
    assert all(producer[k] is True for k in ('stage_metadata_not_long_horizon_projected', 'long_chain_metadata_growth_unqualified',
        'filesystem_allocated_bytes_unmeasured', 'temporary_readback_posterior_products_excluded',
        'array_member_bytes_rehashed', 'runtime_projection_uncalibrated', 'storage_projection_uncalibrated',
        'resource_proxies_are_not_samples'))
    assert all(producer[k] is False for k in ('full_native_likelihood_or_joint_frame_semantics_retested',
        'production_launch_allowed', 'scientific_eligibility', 'posterior_qualified', 'gpu'))
    for p in (args.producer, args.producer_transport, Path(__file__), Path('scripts/ancestral_chain_attempt.py'),
              Path('scripts/reference_measurement_union_sources.py')):
        bind(pins, p)
    # Every derived byte is fresh-checked above; other input bindings remain declared shared dependencies.
    result = dict(status='passed_independent_full_ancestral_derived_storage_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), **checks, resource_proxy_chain_ids=sorted(proxy_rows),
        every_derived_file_rehashed=True, every_member_values_shape_dtype_digest_checked=True,
        every_original_failure_review_and_comparison_retained=True, source_hashes=pins,
        scalar_or_joint_posterior_adequacy_proved=False, scientific_eligibility=False, posterior_qualified=False,
        production_launch_allowed=False, new_native_runs=0, new_predictions=0, gpu=False,
        scope='All derived bytes/NumPy members, attempts, missing-export dispositions and separate resource '
              'arithmetic replayed. Producer streaming NPY bytes and independent NumPy value loads share '
              'NumPy serialization and original closed inputs. Fourteen peer proxies are disk planning only. '
              'No reconstruction of missing draws, full biological revalidation, long-run bound or launch.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
