#!/usr/bin/env python3
"""Measure all current derived ancestral products and explicit long-run scenarios."""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import zipfile

import numpy as np

from ancestral_chain_attempt import sha
from independent_joint_ancestral_frames import ARRAY_NAMES
from reference_measurement_union_sources import bind


SUCCESS = 'full_joint_short_sampler_output_integrity_checked_not_posterior'


def table(path, rows):
    with path.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('output', 'receipt'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists() and not args.receipt.exists()
    pins, expected = {}, {}

    def doc(path, wanted=None):
        path = Path(path)
        digest = sha(path)
        if wanted is not None:
            assert digest == wanted, str(path)
        bind(pins, path, digest)
        return json.loads(path.read_text())

    c6 = doc('metadata/baliphy_scalar_v6_short_sampler_completed_20261004_v1.json')
    assert c6['status'] == 'complete_verified_full_scalar_v6_short_sampler_qualification_v1'
    archive = doc(c6['full_hash_archive'], c6['full_hash_archive_sha256'])
    assert len(archive['source_hashes']) == c6['bound_source_hashes']
    expected.update({str(Path(p).resolve()): h for p, h in archive['source_hashes'].items()})
    for path in ('metadata/baliphy_log_alpha_v10_full_grid_transport_20261006_v1.json',
                 'metadata/baliphy_log_alpha_v10_full_grid_readback_transport_20261006_v1.json',
                 'metadata/current_ancestral_horizon_cost_readback_transport_20261006_v1.json'):
        transport = doc(path)
        assert transport['original_tool_terminal_exit_code'] == 0
        assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
        for p, h in transport['source_hashes'].items():
            key = str(Path(p).resolve())
            assert key not in expected or expected[key] == h
            expected[key] = h
    native_cost = doc('metadata/current_ancestral_horizon_cost_20261006_v1.json',
                      expected[str(Path('metadata/current_ancestral_horizon_cost_20261006_v1.json').resolve())])
    assert native_cost['computational_cost_baselines'] == 1620 and native_cost['proposed_iterations'] == 10000

    def closed(path):
        return doc(path, expected[str(Path(path).resolve())])

    plans = [closed('metadata/baliphy_scalar_v6_sampler_execution_plan_20261004_v1.json'),
             closed('metadata/baliphy_log_alpha_v10_full_grid_plan_20261005_v1.json')]
    sources = []
    for phase, plan in zip(('V6_original', 'V10_comparison'), plans):
        rows = closed(Path(plan['output']) / 'dispositions.json')
        jobs = {j['chain']['chain_id']: j for j in closed(plan['jobs'])}
        assert len(rows) == len(jobs) and {r['chain_id'] for r in rows} == set(jobs)
        sources.append((phase, plan, rows, jobs))
    assert len(sources[0][2]) == 1620 and len(sources[1][2]) == 24
    failed = {r['chain_id'] for r in sources[0][2] if r['exit_code'] != 0}
    assert len(failed) == 24
    assert {j['source_v6_chain_id'] for j in sources[1][3].values()} == failed
    file_rows, member_rows, attempt_rows, derived_leaves = [], [], [], []
    array_paths = set()
    totals = Counter()
    for phase, plan, rows, jobs in sources:
        root = Path(plan['output'])
        assert {p.name for p in (root / 'chains').iterdir()} == {r['chain_id'] + '.json' for r in rows}
        eligible_dirs = {r['chain_id'] for r in rows if r['status'] == SUCCESS}
        assert {p.name for p in (root / 'frames').iterdir()} == eligible_dirs
        for row in sorted(rows, key=lambda r: r['chain_id']):
            cid = row['chain_id']
            assert closed(root / 'chains' / (cid + '.json')) == row
            baseline = phase == 'V10_comparison' or cid not in failed
            original = cid if phase == 'V6_original' else jobs[cid]['source_v6_chain_id']
            assert row['scientific_eligibility'] is row['posterior_qualified'] is False
            frames = row['joint_frames']
            if row['status'] == SUCCESS:
                assert row['exit_code'] == 0 and [f['iteration'] for f in frames] == [0, 10, 20]
                assert {p.name for p in (root / 'frames' / cid).iterdir()} == {'frame-0.npz', 'frame-10.npz', 'frame-20.npz'}
            else:
                assert not frames and not (root / 'frames' / cid).exists()
            role_bytes = raw_bytes = overhead_bytes = maximum_frame = 0
            for frame in frames:
                path = Path(frame['projection_array'])
                assert path == root / 'frames' / cid / ('frame-' + str(frame['iteration']) + '.npz')
                assert path not in array_paths
                array_paths.add(path)
                assert sha(path) == frame['projection_array_sha256'] == expected[str(path.resolve())]
                bind(pins, path, frame['projection_array_sha256'])
                assert set(frame['arrays']) == ARRAY_NAMES
                member_total = data_total = header_total = 0
                with zipfile.ZipFile(path) as saved:
                    infos = saved.infolist()
                    assert len(infos) == len(ARRAY_NAMES) and {i.filename for i in infos} == {k + '.npy' for k in ARRAY_NAMES}
                    for info in infos:
                        key = info.filename[:-4]
                        reference = frame['arrays'][key]
                        assert info.compress_type == zipfile.ZIP_STORED
                        with saved.open(info) as handle:
                            version = np.lib.format.read_magic(handle)
                            assert version == (1, 0)
                            shape, fortran, dtype = np.lib.format.read_array_header_1_0(handle)
                            header_bytes = handle.tell()
                            assert not fortran and not dtype.hasobject
                            assert list(shape) == reference['shape'] and str(dtype) == reference['dtype']
                            count = math.prod(shape) * dtype.itemsize
                            digest, seen = hashlib.sha256(), 0
                            while True:
                                block = handle.read(65536)
                                if not block:
                                    break
                                digest.update(block)
                                seen += len(block)
                            assert seen == count and digest.hexdigest() == reference['sha256']
                            assert info.file_size == header_bytes + count
                        member_rows.append(dict(phase=phase, chain_id=cid, iteration=frame['iteration'],
                            path=str(path), member=key, dtype=str(dtype), shape=json.dumps(list(shape), separators=(',', ':')),
                            raw_array_bytes=count, npy_header_bytes=header_bytes, member_bytes=info.file_size,
                            array_data_sha256=digest.hexdigest(), scientific_eligibility=False, posterior_qualified=False))
                        member_total += info.file_size
                        data_total += count
                        header_total += header_bytes
                size = path.stat().st_size
                zip_overhead = size - member_total
                assert zip_overhead >= 0 and size == data_total + header_total + zip_overhead
                file_rows.append(dict(phase=phase, chain_id=cid, iteration=frame['iteration'], path=str(path),
                    sha256=frame['projection_array_sha256'], bytes=size, raw_array_bytes=data_total,
                    npy_header_bytes=header_total, zip_overhead_bytes=zip_overhead,
                    scientific_eligibility=False, posterior_qualified=False))
                role_bytes += size
                raw_bytes += data_total
                overhead_bytes += header_total + zip_overhead
                maximum_frame = max(maximum_frame, size)
            attempt_rows.append(dict(phase=phase, chain_id=cid, original_v6_chain_id=original,
                family=row['family'], effective_input_group=row['effective_input_group'], prior_label=row['prior_label'],
                chain_role=row['chain_role'], original_integrity_disposition=row['status'], native_exit_code=row['exit_code'],
                baseline_for_cost_scenario=baseline, measured_frames=len(frames), measured_array_file_bytes=role_bytes,
                measured_raw_array_bytes=raw_bytes, measured_serialization_overhead_bytes=overhead_bytes,
                maximum_observed_frame_bytes=maximum_frame,
                projection_basis='own_three_integrity_checked_frames' if frames else 'not_yet_assigned',
                resource_proxy_donor_chain_id='', resource_proxy_frame_bytes=0, projected_array_bytes=0,
                scientifically_imputed=False, scientific_eligibility=False, posterior_qualified=False))
        # Do not descend into immutable native attempts: they have their separate full census.
        leaves = [p for p in root.iterdir() if p.is_file()]
        leaves += [p for p in (root / 'chains').iterdir() if p.is_file()]
        leaves += [p for p in (root / 'frames').glob('*/*') if p.is_file()]
        assert {p.name for p in root.iterdir() if p.is_dir()} == {'attempts', 'chains', 'frames'}
        for path in sorted(leaves):
            digest = sha(path)
            key = str(path.resolve())
            if key in expected:
                assert digest == expected[key]
                custody = 'closed_original_source'
            elif path == Path(c6['full_hash_archive']):
                assert digest == c6['full_hash_archive_sha256']
                custody = 'closed_completion_archive'
            elif phase == 'V10_comparison' and path == root / 'receipt.json':
                custody = 'newly_measured_supplementary_stage_metadata_not_completion_evidence'
            else:
                assert path.suffix == '.lock' and path.stat().st_size == 0
                custody = 'observed_zero_byte_lock_not_completion_evidence'
            bind(pins, path, digest)
            kind = 'projection_array' if path in array_paths else 'role_disposition' if path.parent.name == 'chains' else 'stage_metadata'
            totals[kind] += path.stat().st_size
            derived_leaves.append(dict(phase=phase, kind=kind, path=str(path), sha256=digest,
                bytes=path.stat().st_size, custody=custody))
    baseline_rows = [r for r in attempt_rows if r['baseline_for_cost_scenario']]
    assert len(attempt_rows) == 1644 and len(baseline_rows) == 1620
    assert len({r['original_v6_chain_id'] for r in baseline_rows}) == 1620
    peers = defaultdict(list)
    for row in baseline_rows:
        if row['measured_frames']:
            peers[row['effective_input_group'], row['prior_label']].append(row)
    for row in attempt_rows:
        if not row['baseline_for_cost_scenario']:
            row['projection_basis'] = 'original_native_failure_retained_not_baseline'
        elif row['measured_frames']:
            row['projected_array_bytes'] = (row['measured_array_file_bytes'] * 1001 + 2) // 3
        else:
            donors = peers[row['effective_input_group'], row['prior_label']]
            assert donors, 'Unestimated resource role: no same-input/same-prior donor'
            donor = max(donors, key=lambda r: (r['maximum_observed_frame_bytes'], r['chain_id']))
            row['projection_basis'] = 'same_input_prior_peer_maximum_frame_resource_proxy_only'
            row['resource_proxy_donor_chain_id'] = donor['chain_id']
            row['resource_proxy_frame_bytes'] = donor['maximum_observed_frame_bytes']
            row['projected_array_bytes'] = donor['maximum_observed_frame_bytes'] * 1001
    assert len(file_rows) == 4818 and len(member_rows) == 48180
    assert sum(r['measured_frames'] > 0 for r in baseline_rows) == 1606
    assert sum(not r['measured_frames'] for r in baseline_rows) == 14
    assert len(derived_leaves) == 6477 and Counter(r['kind'] for r in derived_leaves)['projection_array'] == 4818
    args.output.mkdir(parents=True)
    tables = {'attempt_table': ('derived_attempts.tsv', attempt_rows), 'file_table': ('projection_files.tsv', file_rows),
              'member_table': ('projection_members.tsv', member_rows), 'derived_file_table': ('all_derived_files.tsv', derived_leaves)}
    paths = {}
    for label, (name, rows) in tables.items():
        path = args.output / name
        table(path, rows)
        bind(pins, path)
        paths[label] = str(path)
    bind(pins, Path(__file__))
    bind(pins, 'scripts/independent_joint_ancestral_frames.py')
    bind(pins, 'scripts/reference_measurement_union_sources.py')
    bind(pins, 'scripts/ancestral_chain_attempt.py')
    observed_projection = sum(r['projected_array_bytes'] for r in baseline_rows if r['measured_frames'])
    proxy_projection = sum(r['projected_array_bytes'] for r in baseline_rows if not r['measured_frames'])
    native_bytes = native_cost['linear_native_total_with_latent_allowance_bytes']
    result = dict(status='complete_full_ancestral_derived_storage_census_pending_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), measured_attempts=1644, original_roles=1620,
        effective_inputs=135, original_quartets=405, failed_originals_retained=24, comparison_attempts_retained=24,
        computational_baselines=1620, integrity_checked_baselines=1606, review_baselines_without_exports=14,
        measured_projection_files=len(file_rows), measured_array_members=len(member_rows), measured_derived_files=len(derived_leaves),
        measured_bytes_by_kind=dict(totals), measured_derived_bytes=sum(totals.values()),
        measured_raw_array_bytes=sum(r['raw_array_bytes'] for r in file_rows),
        measured_serialization_overhead_bytes=sum(r['npy_header_bytes'] + r['zip_overhead_bytes'] for r in file_rows),
        largest_observed_projection_file_bytes=max(r['bytes'] for r in file_rows),
        proposed_iterations=10000, saved_frame_interval=10, proposed_frames_per_role=1001,
        linear_observed_role_array_projection_bytes=observed_projection,
        resource_only_peer_proxy_projection_bytes=proxy_projection,
        array_projection_with_review_resource_proxies_bytes=observed_projection + proxy_projection,
        native_projection_with_latent_allowance_bytes=native_bytes,
        native_plus_array_projection_with_resource_proxies_bytes=native_bytes + observed_projection + proxy_projection,
        stage_metadata_not_long_horizon_projected=True, long_chain_metadata_growth_unqualified=True,
        filesystem_allocated_bytes_unmeasured=True, temporary_readback_posterior_products_excluded=True,
        full_native_likelihood_or_joint_frame_semantics_retested=False, array_member_bytes_rehashed=True,
        runtime_projection_uncalibrated=True, storage_projection_uncalibrated=True,
        resource_proxies_are_not_samples=True, production_launch_allowed=False,
        scientific_eligibility=False, posterior_qualified=False, new_native_runs=0, gpu=False, new_predictions=0,
        source_hashes=pins, **paths,
        scope='Complete exported projection bytes/members and all current non-native products in both completed '
              'run trees. Every original failure and numeric-review role retained. Derived serialization '
              'is checked against closed arrays; previous biological/frame integrity checks are dependencies. '
              'Fourteen missing export roles receive labelled same-input/same-prior peer resource proxies '
              'only; no missing draw is constructed or admitted. Native-plus-array linear scenario excludes '
              'future metadata growth, temporary products and filesystem overhead; it is not a disk bound '
              'or permission to launch. No new sampler, posterior, GPU or cost.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
