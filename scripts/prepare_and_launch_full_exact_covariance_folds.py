#!/usr/bin/env python3
"""Launch the complete exact dependency/cone proof, retaining counterexamples."""
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

import psutil

from ancestral_chain_attempt import sha
from full_exact_covariance_sources import load
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from reference_measurement_union_sources import verify
from run_full_exact_covariance_folds import STATUS as PRODUCER_STATUS, SUMMARY_FIELDS
from readback_full_exact_covariance_folds import STATUS as READER_STATUS


def main():
    gate_path = Path('metadata/exact_covariance_folds_software_validation_20261003_v1.json')
    execution_path = Path('metadata/exact_covariance_folds_software_execution_20261003_v1.json')
    transport_path = Path('metadata/exact_covariance_folds_software_transport_20261003_v1.json')
    gate, execution, transport = [json.loads(path.read_text()) for path in [gate_path, execution_path, transport_path]]
    assert gate['status'] == 'passed_full_exact_covariance_fold_software_contracts'
    assert (gate['full_real_cohorts'], gate['full_real_logical_cases'], gate['full_real_cohort_row_occurrences']) == (4340, 75188, 34110120)
    assert len(gate['synthetic_scenarios']) == 10 and gate['covariance_and_cone_roundtrips'] == 200
    assert len(gate['malformed_cases_rejected']) == 12 and gate['source_row_serialization_checked']
    assert execution['status'] == 'exited_zero_with_receipt' and execution['exit_code'] == 0
    assert execution['receipt_sha256'] == sha(gate_path)
    assert transport['status'] == 'verified_original_exact_covariance_software_wait_exited_zero'
    assert transport['actual_tool_terminal_exit_code'] == 0 and transport['invocation_id'] == execution['invocation_id']
    for record in [gate, execution, transport]:
        verify(record['source_hashes'])
        if 'artifacts' in record: verify(record['artifacts'])
    assert psutil.virtual_memory().available >= 16 * 2**30 and shutil.disk_usage('.').free >= 64 * 2**30
    output = 'results/phylogeny/full-exact-uniform-covariance-folds-20261003-v1'
    assert not Path(output).exists()
    scope = ('Full exact named covariance-dependency certification for all4340cohorts,75188logical cases '
        'and both loading modes;8680certificates cover68,220,240cohort-row occurrences. '
        'Exact scaled-integer sparse Grams test pair=target+background,gene=pair/2,model=gene. '
        'A failed identity retains its original kernel and integer counterexample. Named folds are conditional '
        'on exact equality,not rank or tolerance. Full nonnegative forward/right-inverse maps preserve the '
        'uniform covariance cone; target/residual and signed-zero/unsigned-four-times-family-intercept '
        'identities are proved for every cohort. CSR multiplication and independent CSC column outer sums '
        'must agree for every case. All13operators/all4340membership files are rehashed against already '
        'closed operator/design archives; broader inherited source files are not freshly rehashed here. '
        'Producer/independent reader/source/artifact/two-original-journal closure required. No case filtering, '
        'variance attribution,nonuniform-weight acceptance,raw/REML model qualification,optimization, '
        'biological acceptance,GPU,new costs or source/attempt changes to original jobs.')
    resources = dict(checked_utc=datetime.now(timezone.utc).isoformat(), cpus=2, memory_gib=16, swap_gib=0,
        blas_threads=1, workers=1, address_space_gib=12, cpu_seconds_per_stage=14400,
        per_file_limit_mib=128, planning_output_gib=4, minimum_free_disk_gib=64,
        expected_certificates=8680, full_source_row_occurrences=68220240,
        global_family_squared_rows=gate['global_family_squared_rows'],
        one_global_family_int64_dense_kernel_gib=gate['global_family_squared_rows'] * 8 / 2**30,
        software_child_cpu_seconds=execution['child_cpu_seconds'], software_child_peak_rss_bytes=execution['child_peak_rss_bytes'],
        available_memory_gib=psutil.virtual_memory().available / 2**30,
        available_disk_gib=shutil.disk_usage('.').free / 2**30, finish_eta=None, gpu=False, new_cost_usd=0,
        scope='Two-CPU/16GiB/no-swap service and12GiB/14400CPU-second/128MiB-file process caps. One worker keeps sparse kernel work bounded; no full dense n-by-n covariance. Global family-squared-row bound describes potential matrix support,not measured peak memory. FourGiB is planning,not an enforced total-byte cap. Runtime is uncalibrated; no finish ETA or automatic retry.')
    own = ['covariance_exact_folds', 'full_exact_covariance_sources', 'run_full_exact_covariance_folds',
        'readback_full_exact_covariance_folds', 'check_covariance_exact_folds', 'run_exact_covariance_software_stage_v1',
        'prepare_and_launch_full_exact_covariance_folds', 'launch_baliphy_reference_sampler_qualification',
        'launch_full_triad_sequence_geometry_followup', 'run_after_verified_dependencies_v2',
        'close_full_triad_sequence_stage', 'record_completed_process_handoffs_v2']
    pins = {f'scripts/{name}.py': sha(f'scripts/{name}.py') for name in own}
    for path in [gate_path, execution_path, transport_path, Path('/usr/bin/prlimit'),
        Path('metadata/exact_covariance_folds_software_resources_20261003_v1.json')]: pins[str(path)] = sha(path)
    pins.update(transport['source_hashes']); pins.update(execution['artifacts'])
    path = 'metadata/full_exact_covariance_folds_plan_20261003_v1.json'
    plan = dict(operator_plan='metadata/full_entity_operator_plan_20261002.json',
        operator_completion='metadata/full_entity_operator_bank_completed_20261002.json',
        design_plan='metadata/full_expanded_model_designs_plan_20261002_v2.json',
        design_completion='metadata/full_expanded_model_designs_v2_completed_20261002.json',
        expected={'cohorts': 4340}, output=output, pins=pins, resources=resources, scope=scope,
        completion='metadata/full_exact_covariance_folds_completed_20261003_v1.json',
        launch_inventory='metadata/full_exact_covariance_folds_launches_20261003_v1.json')
    source_path = 'metadata/full_exact_covariance_folds_prelaunch_source_plan_20261003_v1.json'
    create(source_path, plan); source, bindings = load(plan, Path(source_path)); plan['pins'].update(bindings); create(path, plan)
    prefix = ['/usr/bin/prlimit', '--as=' + str(12 * 2**30), '--cpu=14400', '--fsize=' + str(128 * 2**20), '--']
    producer = launch('full-exact-covariance-folds', prefix + [sys.executable, 'scripts/run_full_exact_covariance_folds.py',
        '--plan', path], [], path, cpus=2, memory=16)
    reader = launch('full-exact-covariance-folds-readback', prefix + [sys.executable, 'scripts/readback_full_exact_covariance_folds.py',
        '--plan', path], [producer], path, cpus=2, memory=16)
    completion_plan = 'metadata/full_exact_covariance_folds_completion_plan_20261003_v1.json'
    create(completion_plan, dict(source_plan=path, producer_receipt=output + '/receipt.json',
        independent_readback=output + '/readback.json', producer_status=PRODUCER_STATUS, reader_status=READER_STATUS,
        completed_status='complete_verified_full_exact_uniform_covariance_folds', summary_fields=SUMMARY_FIELDS,
        launches=[producer, reader], pins={q: sha(q) for q in [path, producer, reader,
            'scripts/close_full_triad_sequence_stage.py', 'scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'], scope=scope))
    closer = launch('full-exact-covariance-folds-closure', [sys.executable, 'scripts/close_full_triad_sequence_stage.py',
        '--plan', completion_plan], [producer, reader], completion_plan, cpus=2, memory=16)
    create(plan['launch_inventory'], dict(status='launched_full_exact_uniform_covariance_fold_proof',
        source_plan=path, source_plan_sha256=sha(path), launches=[producer, reader, closer],
        expected_cohorts=4340, expected_certificates=8680, raw_reml_basis_qualification_complete=False,
        production_fitting_launched=False, all_eight_aims_incomplete=True, gpu=False, new_cost_usd=0))


if __name__ == '__main__': main()
