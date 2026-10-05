#!/usr/bin/env python3
"""Freeze all 24 V10 paired roles after actual job and reader-software closure."""
from datetime import datetime, timezone
import json
from pathlib import Path

import psutil

from baliphy_log_alpha_v10_full_jobs_v1 import sources, validate_jobs
from reference_measurement_union_sources import bind, verify
from run_baliphy_log_alpha_v10_full_grid_v1 import gates


def save(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def main():
    plan_path = Path('metadata/baliphy_log_alpha_v10_full_grid_plan_20261005_v1.json')
    resources_path = Path('metadata/baliphy_log_alpha_v10_full_grid_resources_20261005_v1.json')
    reader_resources_path = Path('metadata/baliphy_log_alpha_v10_full_grid_readback_resources_20261005_v1.json')
    assert not any(p.exists() for p in [plan_path, resources_path, reader_resources_path])
    descriptors = []
    for name, status in [('jobs', 'passed_exact_all24_V10_full_comparison_job_contracts_v1'),
                         ('reader', 'passed_V10_full_reader_closed_trace_and_review_state_integration')]:
        stem = 'metadata/baliphy_log_alpha_v10_full_' + name + '_software'
        descriptors.append(dict(receipt=stem + '_validation_20261005_v1.json',
                                transport=stem + '_transport_20261005_v1.json', status=status))
    pins = gates(dict(gates=descriptors))
    original, source_pins, mapping, original_scope = sources()
    qualified = json.loads(Path(descriptors[0]['receipt']).read_text())
    jobs_path = Path(qualified['jobs'])
    jobs = json.loads(jobs_path.read_text())
    job_scope = validate_jobs(jobs, original, original_scope)
    assert qualified['scope_check'] == job_scope and qualified['mapping'] == mapping
    for path, digest in source_pins.items():
        bind(pins, path, digest)
    old_plan_path = Path('metadata/baliphy_joint_fasta_v7_failure_grid_plan_20261004_v1.json')
    old_plan = json.loads(old_plan_path.read_text())
    old_rows = json.loads((Path(old_plan['output']) / 'dispositions.json').read_text())
    actual_elapsed = [r['elapsed_worker_seconds'] for r in old_rows]
    worker_hours = sum(actual_elapsed) / 3600
    resources = dict(old_plan['resources'])
    resources.update(prepared_utc=datetime.now(timezone.utc).isoformat(),
        available_ram_gib=psutil.virtual_memory().available / 2**30,
        free_disk_gib=psutil.disk_usage('.').free / 2**30,
        completed_V7_native_roles=24, measured_V7_total_worker_hours=worker_hours,
        measured_V7_four_worker_lower_bound_wall_hours=worker_hours / 4,
        measured_V7_role_wall_seconds_range=[min(actual_elapsed), max(actual_elapsed)],
        prospective_wall_hours=[2,12], runtime_uncalibrated_beyond_one_case=False,
        runtime_is_uncalibrated_for_V10=True,
        scope='All24new V10full-input paired roles,4CPUworkers200GiB0swap192GiBnative leases '
              '(4x48GiB)+8GiBcontroller headroom. Original seeds/iterations/input/prior/CPU/wall/file/8MiBstack caps preserved. '
              'Complete previous24role worker time informs planning2-12h, not a measured V10ETA or adequate posterior. '
              'Only native program and new output namespace change; no automatic retry or old chain continuation.')
    assert resources['available_ram_gib'] >= 200 and resources['free_disk_gib'] >= 256
    reader_resources = dict(prepared_utc=resources['prepared_utc'], cpus=2, memory_gib=32, swap_gib=0,
        blas_threads=1, address_space_gib=24, cpu_seconds_per_stage=7200, wall_seconds_per_stage=10800,
        per_file_limit_mib=2048, minimum_available_ram_gib=32, minimum_free_disk_gib=256,
        output_allowance_gib=8, estimated_wall_hours=[0.1,3], runtime_is_uncalibrated=True,
        gpu=False, new_cost_usd=0, new_native_runs=0,
        scope='All24scientific file pairs/latent traces and serialized native/scalar/joint replay only after '
              'actual original full producer API/native/whole-journal/source closure;2CPU32GiB0swap24GiBAS '
              'with8GiBcontroller soft AS. No new native sampler, review-array admission or posterior acceptance.')
    save(resources_path, resources)
    save(reader_resources_path, reader_resources)
    for path in [Path(__file__), jobs_path, Path(mapping), old_plan_path, resources_path, reader_resources_path,
                 *[Path('scripts', name + '.py') for name in [
                     'baliphy_log_alpha_v10_full_jobs_v1', 'baliphy_log_alpha_v10_full_pairs_v1',
                     'run_baliphy_log_alpha_v10_full_grid_v1', 'readback_baliphy_log_alpha_v10_full_grid_v1',
                     'readback_baliphy_log_alpha_v10', 'run_baliphy_scalar_v6_sampler_v1',
                     'reference_sampler_memory_budget']]]:
        bind(pins, path)
    verify(pins)
    root = 'results/ancestral/baliphy-log-alpha-v10-all24-full-input-comparisons-20261005-v1'
    assert not Path(root).exists()
    plan = dict(status='prepared_all24_V10_full_input_comparisons_after_original_software_closure',
        prepared_utc=resources['prepared_utc'], jobs=str(jobs_path.resolve()), mapping=mapping,
        output=root, job_scope=job_scope, gates=descriptors, resources=resources,
        reader_resources=str(reader_resources_path), pins=pins,
        scope='Both full622-tip effective inputs, all3priors/fourchains/24original comparison seeds and unchanged '
              'native limits in fresh attempts; no original failure restarted or discarded. All native outcomes '
              'and every six-file pair retained. A separate qualified strict JSON/90digit Decimal reader checks '
              'every available latent row and serialized native output after complete original producer closure. '
              'Diagnostic errors, native failures, altered scientific outputs and nonfinite parameters remain '
              'explicit reviews. No tolerance relaxation, old reviewed-array admission, historical latent recovery, '
              'adequate ancestral posterior, full1620role V10sampling, GPU prediction or completed biological aim.')
    save(plan_path, plan)
    print(json.dumps({k: v for k, v in plan.items() if k != 'pins'}, indent=2))


if __name__ == '__main__':
    main()
