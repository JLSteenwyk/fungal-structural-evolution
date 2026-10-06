#!/usr/bin/env python3
"""Replay the complete current cost census with separate byte and integer arithmetic."""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind


def physical(path):
    digest = hashlib.sha256()
    size = lines = longest = pending = 0
    with path.open('rb') as handle:
        while True:
            block = handle.read(65536)
            if not block:
                break
            digest.update(block)
            size += len(block)
            pieces = block.split(b'\n')
            if len(pieces) == 1:
                pending += len(block)
                continue
            longest = max(longest, pending + len(pieces[0]) + 1,
                          max((len(p) + 1 for p in pieces[1:-1]), default=0))
            lines += len(pieces) - 1
            pending = len(pieces[-1])
    if pending:
        lines += 1
        longest = max(longest, pending)
    return dict(bytes=size, lines=lines, longest_line_bytes=longest, sha256=digest.hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('producer', 'producer-transport', 'receipt'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    producer = json.loads(args.producer.read_text())
    transport = json.loads(args.producer_transport.read_text())
    assert producer['status'] == 'complete_current_full_ancestral_horizon_cost_census_pending_readback'
    assert transport['original_tool_terminal_exit_code'] == 0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    assert transport['validation_sha256'] == sha(args.producer)
    pins = dict(transport['source_hashes'])
    p6 = json.loads(Path('metadata/baliphy_scalar_v6_sampler_execution_plan_20261004_v1.json').read_text())
    p10 = json.loads(Path('metadata/baliphy_log_alpha_v10_full_grid_plan_20261005_v1.json').read_text())
    d6 = json.loads((Path(p6['output']) / 'dispositions.json').read_text())
    d10 = json.loads((Path(p10['output']) / 'dispositions.json').read_text())
    j6 = {r['chain']['chain_id']: r for r in json.loads(Path(p6['jobs']).read_text())}
    j10 = {r['chain']['chain_id']: r for r in json.loads(Path(p10['jobs']).read_text())}
    source = {('V6_original', r['chain_id']): r for r in d6}
    source.update({('V10_comparison', r['chain_id']): r for r in d10})
    assert len(d6) == len(j6) == 1620 and len(d10) == len(j10) == 24 and len(source) == 1644
    failed = {r['chain_id'] for r in d6 if r['exit_code'] != 0}
    assert len(failed) == 24 and {j['source_v6_chain_id'] for j in j10.values()} == failed
    observations_path = Path('results/ancestral/full-baliphy-scalar-v6-resource-observation-20261004-v1/roles.json')
    observations = {r['chain_id']: r for r in json.loads(observations_path.read_text())}
    assert set(observations) == set(j6)
    with Path(producer['attempt_table']).open() as handle:
        attempts = list(csv.DictReader(handle, delimiter='\t'))
    with Path(producer['file_table']).open() as handle:
        files = list(csv.DictReader(handle, delimiter='\t'))
    by_key = {(r['phase'], r['chain_id']): r for r in attempts}
    assert len(attempts) == len(by_key) == len(source) and set(by_key) == set(source)
    file_groups = defaultdict(dict)
    checked_paths = set()
    byte_totals = Counter()
    maximum_latent_line = projected_bytes = 0
    iterations = producer['proposed_iterations']
    assert iterations == 10000
    for row in files:
        key = row['phase'], row['chain_id']
        assert key in source
        path = Path(row['path'])
        assert str(path) not in checked_paths
        checked_paths.add(str(path))
        measured = physical(path)
        assert measured['sha256'] == row['sha256'] == pins[str(path)]
        for name in ('bytes', 'longest_line_bytes'):
            assert measured[name] == int(row[name])
        assert measured['lines'] == int(row['physical_lines'])
        baseline = row['phase'] == 'V10_comparison' or row['chain_id'] not in failed
        assert row['baseline_for_cost_scenario'] == str(baseline)
        if baseline and path.name in {'C1.P1.fastas', 'C1.P1.site-property-samples.jsonl'}:
            numerator, denominator = 1001, 3
            kind = 'saved_frame_count_linear_scenario'
            if path.name.endswith('.fastas'):
                with path.open('rb') as handle:
                    frame_ids = [int(line.partition(b'=')[2]) for line in handle if line.startswith(b'iterations = ')]
                assert frame_ids == [0, 10, 20]
            else:
                assert measured['lines'] == 3
        elif baseline and path.name in {'C1.log', 'C1.log.json', 'C1.P1.log-alpha-samples.jsonl'}:
            numerator, denominator = 10001, 21
            kind = 'scalar_state_count_linear_scenario'
            assert measured['lines'] == 22
        else:
            numerator, denominator, kind = 1, 1, 'fixed_bytes_scenario'
        quotient, remainder = divmod(measured['bytes'] * numerator, denominator)
        projection = quotient + bool(remainder) if baseline else 0
        assert int(row['projected_bytes']) == projection
        assert int(row['multiplier_numerator']) == numerator and int(row['multiplier_denominator']) == denominator
        assert row['projection_kind'] == kind
        projected_bytes += projection
        byte_totals[path.name] += measured['bytes']
        if path.name == 'C1.P1.log-alpha-samples.jsonl':
            maximum_latent_line = max(maximum_latent_line, measured['longest_line_bytes'])
        assert str(path) not in file_groups[key]
        file_groups[key][str(path)] = row
    baseline_keys = set()
    elapsed = baseline_elapsed = Decimal(0)
    reviews = 0
    for key, expected in sorted(source.items()):
        row = by_key[key]
        phase, cid = key
        job = j6[cid] if phase == 'V6_original' else j10[cid]
        original_id = cid if phase == 'V6_original' else job['source_v6_chain_id']
        baseline = phase == 'V10_comparison' or cid not in failed
        assert row['original_v6_chain_id'] == original_id
        assert row['baseline_for_cost_scenario'] == str(baseline)
        for field in ('family', 'effective_input_group', 'prior_label', 'chain_role', 'native_exit_code', 'source_seed'):
            source_field = {'native_exit_code': 'exit_code', 'source_seed': 'seed'}.get(field, field)
            assert row[field] == str(expected[source_field])
        assert row['original_integrity_disposition'] == expected['status']
        assert row['scientific_eligibility'] == row['posterior_qualified'] == 'False'
        native_path = Path(expected['native_receipt'])
        assert sha(native_path) == expected['native_receipt_sha256']
        native = json.loads(native_path.read_text())
        wanted = {str(native_path.parent / name): h for name, h in native['artifacts'].items()}
        assert set(wanted) == set(file_groups[key])
        assert all(wanted[p] == file_groups[key][p]['sha256'] for p in wanted)
        assert int(row['measured_native_artifact_bytes']) == sum(int(r['bytes']) for r in file_groups[key].values())
        assert int(row['largest_measured_native_file_bytes']) == max(int(r['bytes']) for r in file_groups[key].values())
        assert int(row['projected_existing_file_bytes']) == sum(int(r['projected_bytes']) for r in file_groups[key].values())
        assert int(row['declared_address_space_bytes']) == job['memory_reservation_bytes']
        value = Decimal(str(native['elapsed_seconds']))
        assert Decimal(row['measured_native_worker_wall_seconds']) == value
        elapsed += value
        if baseline:
            assert native['exit_code'] == 0
            baseline_keys.add(original_id)
            baseline_elapsed += value
            reviews += expected['status'] != 'full_joint_short_sampler_output_integrity_checked_not_posterior'
        if phase == 'V6_original':
            obs = observations[original_id]
            assert int(row['observed_v6_vm_hwm_bytes']) == obs['maximum_reported_memory_bytes']['VmHWM']
            assert Decimal(row['observed_v6_cpu_seconds_lower_bound']) == Decimal(str(obs['maximum_observed_cpu_seconds']))
            assert int(row['live_memory_observations']) == obs['live_observations']
        else:
            assert row['observed_v6_vm_hwm_bytes'] == row['observed_v6_cpu_seconds_lower_bound'] == row['live_memory_observations'] == ''
    assert baseline_keys == set(j6) and len(baseline_keys) == 1620
    assert len({r['effective_input_group'] for r in d6}) == 135
    assert len({r['model_input_identity'] for r in d6}) == 405
    assert Counter(r['prior_label'] for r in d6) == {'broad': 540, 'centered': 540, 'package': 540}
    for field, value in dict(original_roles=1620, effective_inputs=135, original_quartets=405,
                             original_native_failures=24, separate_comparison_attempts=24,
                             measured_native_attempts=1644, computational_cost_baselines=1620,
                             proposed_saved_frame_interval=10, expected_scalar_states_per_role=10001,
                             expected_saved_frames_per_role=1001).items():
        assert producer[field] == value
    assert len(files) == producer['measured_native_files'] and dict(byte_totals) == producer['measured_bytes_by_filename']
    assert sum(byte_totals.values()) == producer['measured_native_artifact_bytes']
    assert elapsed == Decimal(producer['measured_consumed_native_worker_wall_seconds'])
    assert baseline_elapsed == Decimal(producer['baseline_native_worker_wall_seconds'])
    assert projected_bytes == producer['linear_projected_existing_native_bytes']
    assert maximum_latent_line == producer['largest_observed_latent_physical_line_bytes']
    extra = 1596 * 10002 * maximum_latent_line
    assert extra == producer['missing_latent_output_allowance_bytes']
    assert projected_bytes + extra == producer['linear_native_total_with_latent_allowance_bytes']
    assert reviews == producer['scalar_review_roles_in_cost_baselines']
    with localcontext() as context:
        context.prec = 50
        for row, workers in zip(producer['worker_scenarios'], [4, 8, 16]):
            assert row['workers'] == workers
            value = baseline_elapsed * 500
            assert Decimal(row['native_worker_wall_hours']) == value / 3600
            assert Decimal(row['ideal_worker_only_wall_days']) == value / (86400 * workers)
        historical = json.loads(Path('metadata/baliphy_horizon_resources_completed_20261003.json').read_text())
        assert Decimal(producer['historical_1000_iteration_successful_only_projected_worker_hours']) == \
            Decimal(str(historical['linear_successful_worker_seconds'])) / 3600
        assert producer['historical_failed_roles_unestimated'] == historical['selected_failed_chains']
    telemetry = json.loads(Path('metadata/baliphy_scalar_v6_resource_observation_completed_20261004_v1.json').read_text())
    assert producer['original_v6_maximum_observed_cgroup_bytes'] == telemetry['maximum_observed_cgroup_memory_bytes']
    assert producer['original_v6_maximum_reported_cgroup_peak_bytes'] == telemetry['maximum_observed_cgroup_reported_peak_bytes']
    assert producer['production_launch_allowed'] is producer['scientific_eligibility'] is producer['posterior_qualified'] is False
    assert producer['runtime_projection_uncalibrated'] and producer['storage_projection_uncalibrated']
    assert producer['v10_per_native_live_memory_unobserved'] and producer['long_horizon_peak_memory_unqualified']
    # Every declared source/output binding closes, including tables and compact source proofs.
    for path, digest in pins.items():
        if path not in checked_paths:
            assert sha(path) == digest, path
    for path in (args.producer, args.producer_transport, Path(__file__)):
        bind(pins, path)
    result = dict(status='passed_independent_complete_current_ancestral_horizon_cost_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), original_roles=1620,
        native_attempts_checked=1644, native_files_checked=len(files), native_bytes_checked=sum(byte_totals.values()),
        cost_baselines_checked=1620, original_failed_roles_retained=24, separate_comparisons_retained=24,
        native_file_measurements_reconstructed=True, source_attempt_dispositions_reconstructed=True,
        exact_integer_storage_scenarios_reconstructed=True, decimal_worker_scenarios_reconstructed=True,
        producer_sha256=sha(args.producer), source_hashes=pins,
        production_launch_allowed=False, runtime_projection_uncalibrated=True, storage_projection_uncalibrated=True,
        posterior_qualified=False, scientific_eligibility=False, new_native_runs=0, gpu=False,
        scope='Separate chunked byte/line counters, divmod ceiling arithmetic and Decimal sums replay all 1644 '
              'attempts and every native file, original failures, comparison aliases, memory observation sources '
              'and full 1620-role cost scenarios. Underlying native sources/closure/SHA helpers are shared. '
              'Reproducing a scenario does not establish linear scaling, future memory safety, convergence or biology.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
