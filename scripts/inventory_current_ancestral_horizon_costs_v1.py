#!/usr/bin/env python3
"""Measure every current ancestral native attempt before longer-horizon planning."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind


SAMPLED = {'C1.P1.fastas', 'C1.P1.site-property-samples.jsonl'}
SCALAR = {'C1.log', 'C1.log.json', 'C1.P1.log-alpha-samples.jsonl'}


def measure(path):
    digest = hashlib.sha256()
    size = lines = longest = 0
    frame_iterations = []
    with path.open('rb') as handle:
        for line in handle:
            digest.update(line)
            size += len(line)
            lines += 1
            longest = max(longest, len(line))
            if line.startswith(b'iterations = '):
                frame_iterations.append(int(line.split(b'=', 1)[1]))
    return dict(bytes=size, sha256=digest.hexdigest(), lines=lines,
                longest_line_bytes=longest, fasta_frame_iterations=frame_iterations)


def write_table(path, rows):
    with path.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--iterations', type=int, default=10000)
    args = parser.parse_args()
    assert args.iterations > 1000 and args.iterations % 10 == 0
    assert not args.output.exists() and not args.receipt.exists()
    pins, expected = {}, {}

    def doc(path, digest=None):
        path = Path(path)
        actual = sha(path)
        if digest is not None:
            assert actual == digest, str(path)
        bind(pins, path, actual)
        return json.loads(path.read_text())

    def archive(completion):
        proof = doc(completion['full_hash_archive'], completion['full_hash_archive_sha256'])
        assert len(proof['source_hashes']) == completion['bound_source_hashes']
        for p, h in proof['source_hashes'].items():
            key = str(Path(p).resolve())
            assert key not in expected or expected[key] == h
            expected[key] = h

    def closed(path, digest=None):
        key = str(Path(path).resolve())
        assert key in expected, key
        if digest is not None:
            assert expected[key] == digest
        return doc(path, expected[key])

    c6 = doc('metadata/baliphy_scalar_v6_short_sampler_completed_20261004_v1.json')
    assert c6['status'] == 'complete_verified_full_scalar_v6_short_sampler_qualification_v1'
    archive(c6)
    telemetry = doc('metadata/baliphy_scalar_v6_resource_observation_completed_20261004_v1.json')
    assert telemetry['status'] == 'complete_verified_full_scalar_v6_sampler_resource_observation_v1'
    archive(telemetry)
    p6 = closed('metadata/baliphy_scalar_v6_sampler_execution_plan_20261004_v1.json')
    jobs6 = closed(p6['jobs'])
    root6 = Path(p6['output'])
    dispositions6 = closed(root6 / 'dispositions.json')
    producer6 = closed(c6['producer_receipt'], c6['producer_receipt_sha256'])
    assert producer6['full_chains'] == len(jobs6) == len(dispositions6) == 1620
    telemetry_producer = closed(telemetry['producer_receipt'], telemetry['producer_receipt_sha256'])
    role_observations_path = Path(telemetry['producer_receipt']).parent / 'roles.json'
    observations = closed(role_observations_path)
    observations = {r['chain_id']: r for r in observations}
    jobs6 = {j['chain']['chain_id']: j for j in jobs6}
    assert set(jobs6) == set(observations) == {r['chain_id'] for r in dispositions6}

    p10 = doc('metadata/baliphy_log_alpha_v10_full_grid_plan_20261005_v1.json')
    r10 = doc('metadata/baliphy_log_alpha_v10_full_grid_readback_20261006_v1.json')
    t10 = doc('metadata/baliphy_log_alpha_v10_full_grid_readback_transport_20261006_v1.json')
    assert t10['original_tool_terminal_exit_code'] == 0
    assert t10['validation_sha256'] == sha('metadata/baliphy_log_alpha_v10_full_grid_readback_20261006_v1.json')
    assert r10['full_input_logger_noninterference_passed'] and r10['full_latent_trace_integrity_passed']
    assert r10['roles'] == 24 and r10['diagnostic_rows'] == 504
    for p, h in t10['source_hashes'].items():
        key = str(Path(p).resolve())
        assert key not in expected or expected[key] == h
        expected[key] = h
    jobs10 = closed(p10['jobs'])
    rows10 = closed(Path(p10['output']) / 'dispositions.json')
    jobs10 = {j['chain']['chain_id']: j for j in jobs10}
    comparisons = {j['source_v6_chain_id']: j['chain']['chain_id'] for j in jobs10.values()}
    failed = {r['chain_id'] for r in dispositions6 if r['exit_code'] != 0}
    assert len(failed) == len(comparisons) == 24 and set(comparisons) == failed
    assert all(r['exit_code'] == 0 for r in rows10)
    old_horizon = doc('metadata/baliphy_horizon_resources_completed_20261003.json')
    assert old_horizon['full_chains'] == 1620 and old_horizon['proposed_iterations'] == args.iterations

    files, attempts = [], []
    byte_totals = Counter()
    latent_maximum_line = 0
    for phase, dispositions, jobs in [('V6_original', dispositions6, jobs6), ('V10_comparison', rows10, jobs10)]:
        for row in sorted(dispositions, key=lambda r: r['chain_id']):
            cid = row['chain_id']
            job = jobs[cid]
            original_id = cid if phase == 'V6_original' else job['source_v6_chain_id']
            baseline = (phase == 'V6_original' and cid not in failed) or phase == 'V10_comparison'
            native_path = Path(row['native_receipt'])
            native = closed(native_path, row['native_receipt_sha256'])
            assert native['exit_code'] == row['exit_code']
            assert Decimal(str(native['elapsed_seconds'])) == Decimal(str(row['elapsed_worker_seconds']))
            measured_bytes = projected_bytes = largest = 0
            counters = Counter()
            for name, digest in sorted(native['artifacts'].items()):
                path = native_path.parent / name
                assert expected[str(path.resolve())] == digest
                measured = measure(path)
                assert measured['sha256'] == digest and measured['bytes'] == path.stat().st_size
                bind(pins, path, digest)
                filename = path.name
                measured_bytes += measured['bytes']
                largest = max(largest, measured['bytes'])
                byte_totals[filename] += measured['bytes']
                if baseline and filename in SAMPLED:
                    if filename == 'C1.P1.fastas':
                        assert measured['fasta_frame_iterations'] == [0, 10, 20]
                    else:
                        assert measured['lines'] == 3
                    numerator, denominator, kind = args.iterations // 10 + 1, 3, 'saved_frame_count_linear_scenario'
                elif baseline and filename in SCALAR:
                    assert measured['lines'] == 22
                    numerator, denominator, kind = args.iterations + 1, 21, 'scalar_state_count_linear_scenario'
                else:
                    numerator, denominator, kind = 1, 1, 'fixed_bytes_scenario'
                projection = (measured['bytes'] * numerator + denominator - 1) // denominator if baseline else 0
                projected_bytes += projection
                counters[kind] += 1
                if filename == 'C1.P1.log-alpha-samples.jsonl':
                    latent_maximum_line = max(latent_maximum_line, measured['longest_line_bytes'])
                files.append(dict(phase=phase, chain_id=cid, original_v6_chain_id=original_id,
                    baseline_for_cost_scenario=baseline, native_exit_code=row['exit_code'],
                    path=str(path), sha256=digest, bytes=measured['bytes'], physical_lines=measured['lines'],
                    longest_line_bytes=measured['longest_line_bytes'], projection_kind=kind,
                    multiplier_numerator=numerator, multiplier_denominator=denominator,
                    projected_bytes=projection))
            observed = observations[original_id] if phase == 'V6_original' else None
            attempts.append(dict(phase=phase, chain_id=cid, original_v6_chain_id=original_id,
                family=row['family'], effective_input_group=row['effective_input_group'],
                prior_label=row['prior_label'], chain_role=row['chain_role'], source_seed=row['seed'],
                native_exit_code=row['exit_code'], original_integrity_disposition=row['status'],
                baseline_for_cost_scenario=baseline, measured_native_worker_wall_seconds=str(native['elapsed_seconds']),
                measured_native_artifact_bytes=measured_bytes, largest_measured_native_file_bytes=largest,
                projected_existing_file_bytes=projected_bytes, declared_address_space_bytes=job['memory_reservation_bytes'],
                observed_v6_vm_hwm_bytes=observed['maximum_reported_memory_bytes']['VmHWM'] if observed else '',
                observed_v6_cpu_seconds_lower_bound=str(observed['maximum_observed_cpu_seconds']) if observed else '',
                live_memory_observations=observed['live_observations'] if observed else '',
                scientific_eligibility=False, posterior_qualified=False))
    baseline_rows = [r for r in attempts if r['baseline_for_cost_scenario']]
    assert len(attempts) == 1644 and len(baseline_rows) == 1620
    assert len({r['original_v6_chain_id'] for r in baseline_rows}) == 1620
    assert len({r['effective_input_group'] for r in baseline_rows}) == 135
    assert Counter(r['prior_label'] for r in baseline_rows) == {'broad': 540, 'centered': 540, 'package': 540}
    assert all(r['native_exit_code'] == 0 for r in baseline_rows)
    baseline_wall = sum((Decimal(r['measured_native_worker_wall_seconds']) for r in baseline_rows), Decimal(0))
    consumed_wall = sum((Decimal(r['measured_native_worker_wall_seconds']) for r in attempts), Decimal(0))
    existing_projection = sum(r['projected_existing_file_bytes'] for r in baseline_rows)
    latent_allowance = 1596 * (args.iterations + 2) * latent_maximum_line
    with localcontext() as context:
        context.prec = 50
        projected_wall = baseline_wall * args.iterations / 20
        scenarios = [dict(workers=n, native_worker_wall_hours=str(projected_wall / 3600),
                          ideal_worker_only_wall_days=str(projected_wall / (3600 * 24 * n))) for n in [4, 8, 16]]
        historical_hours = str(Decimal(str(old_horizon['linear_successful_worker_seconds'])) / 3600)
    args.output.mkdir(parents=True)
    file_table = args.output / 'native_files.tsv'
    attempt_table = args.output / 'native_attempts.tsv'
    write_table(file_table, files)
    write_table(attempt_table, attempts)
    for p in [Path(__file__), file_table, attempt_table]:
        bind(pins, p)
    result = dict(status='complete_current_full_ancestral_horizon_cost_census_pending_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), original_roles=1620, effective_inputs=135,
        original_quartets=405, original_native_failures=24, separate_comparison_attempts=24,
        measured_native_attempts=1644, computational_cost_baselines=1620, measured_native_files=len(files),
        original_native_failures_retained=True, comparisons_are_not_biological_replicates=True,
        measured_bytes_by_filename=dict(byte_totals), measured_native_artifact_bytes=sum(byte_totals.values()),
        measured_consumed_native_worker_wall_seconds=str(consumed_wall),
        baseline_native_worker_wall_seconds=str(baseline_wall),
        proposed_iterations=args.iterations, proposed_saved_frame_interval=10,
        expected_scalar_states_per_role=args.iterations+1, expected_saved_frames_per_role=args.iterations//10+1,
        linear_projected_existing_native_bytes=existing_projection,
        missing_latent_output_allowance_bytes=latent_allowance,
        largest_observed_latent_physical_line_bytes=latent_maximum_line,
        linear_native_total_with_latent_allowance_bytes=existing_projection+latent_allowance,
        worker_scenarios=scenarios, historical_1000_iteration_successful_only_projected_worker_hours=historical_hours,
        historical_failed_roles_unestimated=old_horizon['selected_failed_chains'],
        original_v6_maximum_observed_cgroup_bytes=telemetry['maximum_observed_cgroup_memory_bytes'],
        original_v6_maximum_reported_cgroup_peak_bytes=telemetry['maximum_observed_cgroup_reported_peak_bytes'],
        v10_per_native_live_memory_unobserved=True, long_horizon_peak_memory_unqualified=True,
        scalar_review_roles_in_cost_baselines=sum(r['original_integrity_disposition'] !=
            'full_joint_short_sampler_output_integrity_checked_not_posterior' for r in baseline_rows),
        file_table=str(file_table), attempt_table=str(attempt_table), source_hashes=pins,
        runtime_projection_uncalibrated=True, storage_projection_uncalibrated=True,
        production_launch_allowed=False, posterior_qualified=False, scientific_eligibility=False,
        goal_finish_eta=None, convergence_eta=None, new_native_runs=0, gpu=False, new_predictions=0,
        scope='Complete current 1620-role computational cost census with all 24 original failures and 24 separate '
              'V10 comparisons retained. Every used native artifact rehashed and sized; observed memory remains '
              'sampling telemetry, not a guaranteed peak. Linear frame/state-count byte and worker-wall scenarios '
              'are planning assumptions, not bounds, new validated biological samples, convergence or an ETA. '
              'Historical 1000-iteration costs are explicitly older-code/successful-only reference evidence. '
              'Future full long-horizon construction, parsers, calibration and safe caps remain required.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
