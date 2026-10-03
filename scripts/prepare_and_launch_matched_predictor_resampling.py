#!/usr/bin/env python3
"""Launch all paired site/block draws after closed inputs and software proof."""
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

import psutil

from ancestral_chain_attempt import sha
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from matched_predictor_branch_inputs import verify
from matched_predictor_resampling import load, add_splits, MODES, SUMMARY_FIELDS
from run_matched_predictor_resampling import STATUS as PRODUCER_STATUS
from readback_matched_predictor_resampling import STATUS as READER_STATUS


def main():
    validation = Path('metadata/matched_predictor_resampling_software_validation_20261003_v1.json')
    execution_path = Path('metadata/matched_predictor_resampling_software_execution_20261003_v1.json')
    transport_path = Path('metadata/matched_predictor_resampling_software_transport_20261003_v1.json')
    gate = json.loads(validation.read_text())
    assert gate['status'] == 'passed_full_matched_predictor_paired_resampling_software_contracts'
    assert (gate['full_real_inputs'], gate['full_real_resampling_cases'], gate['future_native_roles']) == (133, 53200, 372400)
    assert gate['full_draw_serialization_checked'] and gate['replicates_per_mode'] == 200 and gate['modes'] == MODES
    assert gate['synthetic_native_roles'] == 21 and gate['synthetic_unsuccessful_roles_retained'] == 7
    assert gate['artificial_failed_and_review_array_dispositions'] == 2 and len(gate['malformed_cases_rejected']) == 16
    execution = json.loads(execution_path.read_text()); transport = json.loads(transport_path.read_text())
    assert execution['status'] == 'exited_zero_with_receipt' and execution['exit_code'] == 0
    assert execution['receipt_sha256'] == sha(validation)
    assert transport['status'] == 'verified_original_resampling_software_wait_exited_zero'
    assert transport['actual_tool_terminal_exit_code'] == 0 and transport['invocation_id'] == execution['invocation_id']
    for record in [gate, execution, transport]:
        verify(record['source_hashes'])
        if 'artifacts' in record: verify(record['artifacts'])
    point_plan_path = Path('metadata/matched_predictor_branch_fits_plan_20261003_v1.json')
    point_plan = json.loads(point_plan_path.read_text())
    assert psutil.virtual_memory().available >= 24 * 2**30 and shutil.disk_usage('.').free >= 256 * 2**30
    output = 'results/phylogeny/matched-predictor-paired-resampling-20261003-v1'
    assert not Path(output).exists()
    point_root = Path(point_plan['output']); point_receipt = json.loads((point_root / 'receipt.json').read_text())
    point_files = [p for p in point_root.rglob('*') if p.is_file()]
    point_bytes = sum(p.stat().st_size for p in point_files)
    scope = ('Full paired predictor uncertainty control: every 133 ready fixed-topology input from the closed '
        '8750-case/125-marker/70-view grid, with all 3780 insufficient comparisons preserved upstream. '
        'There are 71 markers and 21 fungal taxa in selected predictor overlap, not full fungal diversity. '
        'Two hundred iid-site and two hundred circular-ten-retained-column-block draws per input, using '
        'identical column indices for AA, AlphaFold and ESMFold states and all taxa; predictions, original '
        'taxon/missingness eligibility, tree topology and state models remain fixed. The 71 identical '
        'alignment/taxon/column draw groups share draws across alternative topologies. All 53200 cases '
        'and 372400 native roles are retained, including failed/invalid cases; no retry, new eligibility '
        'filter or successful-only interval acceptance. Independent full archive/tree/report/array readback '
        'and original producer/reader journal closure required. Branch lengths are state substitutions/site, '
        'not physical displacement or rates/year. Blocks use retained-column order and circular wrapping, '
        'not physical residue distance or guaranteed independent blocks. These are conditional uncertainty '
        'and sensitivity analyses, not predictive-source/model/orthology/framework acceptance, multiple-test '
        'calibration or any completed biological aim. No GPU prediction, costs or old job/source changes.')
    resources = dict(checked_utc=datetime.now(timezone.utc).isoformat(), cpus=8, workers=8, memory_gib=24,
        swap_gib=0, native_threads=1, blas_threads=1, native_address_space_gib=2,
        native_cpu_seconds=300, native_wall_seconds=600, native_per_file_limit_mib=8,
        maximum_simultaneous_native_roles=8, maximum_native_address_space_reservations_gib=16,
        planning_output_gib=128, minimum_free_disk_gib=256, planning_retained_case_files=159600,
        point_fit_total_files=len(point_files), point_fit_total_bytes=point_bytes,
        point_fit_native_seconds_sum=point_receipt['native_seconds_sum'],
        point_scaled_raw_output_gib=point_bytes * 400 / 2**30,
        point_scaled_native_worker_wall_hours_at_eight_workers=point_receipt['native_seconds_sum'] * 400 / 8 / 3600,
        software_child_cpu_seconds=execution['child_cpu_seconds'], software_child_peak_rss_bytes=execution['child_peak_rss_bytes'],
        available_memory_gib=psutil.virtual_memory().available / 2**30,
        available_disk_gib=shutil.disk_usage('.').free / 2**30,
        available_inodes=__import__('os').statvfs('.').f_favail, finish_eta=None, gpu=False, new_cost_usd=0,
        scope='Eight-worker native concurrency and24GiB/no-swap cgroup limits are enforced;2GiB/300CPU-second/600wall-second/8MiB-file native limits terminate and retain failures.128GiB is an output planning allowance, not a hard total-output cap.256GiB free-space guard is checked at start and completed cases. Point-worker scaling excludes input/packing/controller/readback overhead and changed optimization behavior; it is not a finish ETA. Archives reduce raw native file inode count while retaining all bytes.')
    own = ['matched_predictor_resampling', 'run_matched_predictor_resampling', 'readback_matched_predictor_resampling',
        'check_matched_predictor_resampling', 'run_matched_predictor_resampling_stage_v1',
        'prepare_and_launch_matched_predictor_resampling', 'matched_predictor_branch_fits',
        'readback_matched_predictor_branch_fits', 'ancestral_chain_attempt',
        'launch_baliphy_reference_sampler_qualification', 'launch_full_triad_sequence_geometry_followup',
        'close_full_triad_sequence_stage', 'record_completed_process_handoffs_v2', 'run_after_verified_dependencies_v2']
    pins = {f'scripts/{name}.py': sha(f'scripts/{name}.py') for name in own}
    for path in [validation, execution_path, transport_path, point_plan_path,
        Path('metadata/matched_predictor_resampling_software_resources_20261003_v1.json'),
        Path(point_plan['executable']), *map(Path, point_plan['models'].values())]: pins[str(path)] = sha(path)
    pins.update(transport['source_hashes']); pins.update(execution['artifacts'])
    path = 'metadata/matched_predictor_resampling_plan_20261003_v1.json'
    plan = dict(input_plan=point_plan['input_plan'], input_completion=point_plan['input_completion'],
        point_completion='metadata/matched_predictor_branch_fits_completed_20261003_v1.json', point_output=str(point_root),
        executable=point_plan['executable'], models=point_plan['models'], pins=pins, output=output,
        replicates_per_mode=200, modes=MODES, resources=resources, scope=scope,
        completion='metadata/matched_predictor_resampling_completed_20261003_v1.json',
        launch_inventory='metadata/matched_predictor_resampling_launches_20261003_v1.json')
    source_path = 'metadata/matched_predictor_resampling_prelaunch_source_plan_20261003_v1.json'
    create(source_path, plan)
    source, bindings = load(plan, Path(source_path)); add_splits(source, point_root)
    assert len(source['configs']) == 133 and len(set(source['resampling_groups'].values())) == 71
    plan['pins'].update(bindings); create(path, plan)
    producer = launch('matched-predictor-resampling', [sys.executable, 'scripts/run_matched_predictor_resampling.py',
        '--plan', path], [], path, cpus=8, memory=24)
    reader = launch('matched-predictor-resampling-readback', [sys.executable, 'scripts/readback_matched_predictor_resampling.py',
        '--plan', path], [producer], path, cpus=2, memory=16)
    completion_plan = 'metadata/matched_predictor_resampling_completion_plan_20261003_v1.json'
    spec = dict(source_plan=path, producer_receipt=output + '/receipt.json', independent_readback=output + '/readback.json',
        producer_status=PRODUCER_STATUS, reader_status=READER_STATUS,
        completed_status='complete_verified_full_matched_predictor_paired_resampling', summary_fields=SUMMARY_FIELDS,
        launches=[producer, reader], pins={q: sha(q) for q in [path, producer, reader,
        'scripts/close_full_triad_sequence_stage.py', 'scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'], scope=scope)
    create(completion_plan, spec)
    closer = launch('matched-predictor-resampling-closure', [sys.executable, 'scripts/close_full_triad_sequence_stage.py',
        '--plan', completion_plan], [producer, reader], completion_plan, cpus=2, memory=16)
    create(plan['launch_inventory'], dict(status='launched_full_matched_predictor_paired_resampling', source_plan=path,
        source_plan_sha256=sha(path), launches=[producer, reader, closer], expected_cases=53200, expected_native_roles=372400,
        all_eight_aims_incomplete=True, gpu=False, new_cost_usd=0))


if __name__ == '__main__':
    main()
