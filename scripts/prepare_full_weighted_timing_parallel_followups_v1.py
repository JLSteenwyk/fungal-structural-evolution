#!/usr/bin/env python3
"""Bind full parallel timing reader and closure to their original live launches."""
import argparse
import json
from pathlib import Path
import sys

from ancestral_chain_attempt import sha
from full_weighted_timing_parallel_v1 import PRODUCER, READER, SUMMARY
from reference_measurement_union_sources import verify


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--role', choices=['reader', 'closure', 'inventory'], required=True)
    args = parser.parse_args()
    source_plan = Path('metadata/full_weighted_timing_parallel_plan_20261004_v1.json')
    plan = json.loads(source_plan.read_text())
    verify(plan['pins'])
    producer = 'metadata/full_weighted_timing_parallel_producer_launch_20261004_v1.json'
    reader = 'metadata/full_weighted_timing_parallel_reader_launch_20261004_v1.json'
    closer = 'metadata/full_weighted_timing_parallel_closure_launch_20261004_v1.json'
    if args.role == 'inventory':
        launches = [producer, reader, closer]
        records = [json.loads(Path(path).read_text()) for path in launches]
        result = dict(status='queued_complete_full_parallel_timing_after_original_numerical_closure',
            source_plan=str(source_plan), source_plan_sha256=sha(source_plan), launches=launches,
            original_tool_sessions=[r['original_tool_session_id'] for r in records],
            candidate_rows=20832000, maximum_eligible_groups=694400, variance_points_per_group=2,
            independent_numeric_replay='every_selected_group',
            numerical_dependency=plan['numerical_closure_launch'],
            scientific_eligibility=False, weighted_working_model_fits_launched=False,
            gpu=False, original_serial_jobs_restarted=False, all_eight_aims_incomplete=True)
        write('metadata/full_weighted_timing_parallel_launches_20261004_v1.json', result)
        print(json.dumps(result, indent=2))
        return
    if args.role == 'reader':
        dependencies = [producer]
        source = source_plan
        command = ['/usr/bin/prlimit', '--as=' + str(8*2**30) + ':' + str(24*2**30),
            '--cpu=4838400', '--fsize=' + str(2*2**30), '--', sys.executable,
            'scripts/full_weighted_timing_parallel_v1.py', '--plan', str(source_plan), '--reader']
    else:
        dependencies = [producer, reader]
        source = Path('metadata/full_weighted_timing_parallel_completion_plan_20261004_v1.json')
        root = Path(plan['output'])
        write(source, dict(source_plan=str(source_plan), producer_receipt=str(root/'receipt.json'),
            independent_readback=str(root/'readback.json'), producer_status=PRODUCER,
            reader_status=READER, completed_status='complete_verified_full_four_control_input_timing_accounting_v2',
            summary_fields=SUMMARY + ['parallel_workers', 'worker_address_space_gib',
                                     'cached_numeric_inputs_preserved_for_all_cohorts'],
            launches=dependencies, output=plan['completion'], scope=plan['scope'],
            pins={str(path): sha(path) for path in [source_plan, *map(Path, dependencies),
                Path('scripts/close_full_weighted_timing_parallel_v1.py'),
                Path('scripts/record_completed_process_handoffs_v2.py'), Path(__file__)]}))
        command = ['/usr/bin/prlimit', '--as=' + str(24*2**30), '--cpu=7200',
            '--fsize=' + str(2*2**30), '--', sys.executable,
            'scripts/close_full_weighted_timing_parallel_v1.py', '--plan', str(source)]
    wait_plan = Path('metadata/full_weighted_timing_parallel_' + args.role + '_wait_plan_20261004_v1.json')
    write(wait_plan, dict(dependencies=dependencies, command=command,
        pins={str(path): sha(path) for path in [source, *map(Path, dependencies),
            Path('scripts/run_after_verified_dependencies_v2.py'), Path(__file__)]}))
    print(json.dumps(dict(role=args.role, wait_plan=str(wait_plan), command=command), indent=2))


if __name__ == '__main__':
    main()
