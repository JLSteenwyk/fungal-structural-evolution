#!/usr/bin/env python3
"""Queue the full separate scalar/length numerical check and its accounting closure."""
import argparse
import json
from pathlib import Path
import shutil
import sys

import psutil
from independent_baliphy_scalar_sources import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import launch, create
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args(); plan = json.loads(args.plan.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    for key, expected in [('numerical_validation', 'passed_independent_ancestral_scalar_direct_lag_contracts'),
                          ('grid_validation', 'passed_full_independent_baliphy_scalar_grid_and_checkpoint_contracts')]:
        record = json.loads(Path(plan[key]).read_text()); assert record['status'] == expected
        for path, digest in record['source_hashes'].items():
            assert sha(path) == digest
    assert psutil.virtual_memory().available >= 32 * 2**30
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    root = Path(plan['output']); assert not root.exists()
    dependency = json.loads(Path(plan['diagnostic_closure_launch']).read_text())
    dependency['launch'] = plan['diagnostic_closure_launch']
    assert fingerprint(dependency) is None; journal_terminal(dependency)
    producer = launch('independent-baliphy-scalar', 'metadata/independent_baliphy_scalar_wait_plan_20261002.json',
        [sys.executable, 'scripts/prepare_independent_baliphy_scalar_readback.py', '--plan', str(args.plan)],
        [plan['diagnostic_closure_launch']], str(args.plan))
    reader = launch('independent-baliphy-scalar-readback', 'metadata/independent_baliphy_scalar_readback_wait_plan_20261002.json',
        [sys.executable, 'scripts/readback_independent_baliphy_scalar_readback.py', '--plan', str(args.plan), '--output', str(root / 'readback.json')],
        [producer], str(args.plan))
    close_path = plan['completion_plan']
    create(close_path, dict(source_plan=str(args.plan), producer_receipt=str(root / 'receipt.json'),
        independent_readback=str(root / 'readback.json'),
        producer_status='complete_full_independent_baliphy_scalar_comparison_pending_readback',
        reader_status='passed_full_independent_baliphy_scalar_serialized_comparison_readback',
        completed_status='complete_verified_full_independent_baliphy_scalar_comparison',
        summary_fields=SUMMARY_FIELDS, launches=[producer, reader],
        pins={p: sha(p) for p in [str(args.plan), producer, reader,
            'scripts/close_full_triad_sequence_stage.py', 'scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'], scope=plan['scope']))
    closure = launch('independent-baliphy-scalar-closure', 'metadata/independent_baliphy_scalar_closure_wait_plan_20261002.json',
        [sys.executable, 'scripts/close_full_triad_sequence_stage.py', '--plan', close_path], [producer, reader], close_path)
    create(plan['launch_inventory'], dict(status='queued_full_independent_baliphy_scalar_numerical_comparison',
        source_plan=str(args.plan), source_plan_sha256=sha(args.plan), launches=[producer, reader, closure],
        native_sampling_restarted=False, gpu=False, new_cost_usd=0, scope=plan['scope']))


if __name__ == '__main__':
    main()
