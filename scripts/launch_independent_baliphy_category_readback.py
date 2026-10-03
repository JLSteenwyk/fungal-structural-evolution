#!/usr/bin/env python3
"""Launch full categorical replay/readback/closure after pinned software checks."""
import argparse
import json
from pathlib import Path
import shutil
import sys

import psutil

from independent_baliphy_category_sources import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import launch, create
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args(); plan = json.loads(args.plan.read_text())
    assert plan['comparison_tolerances'] == dict(atol=1e-8, rtol=1e-8)
    for path, digest in plan['pins'].items(): assert sha(path) == digest, path
    for key, status in [('numerical_validation', 'passed_independent_ancestral_categorical_contracts'),
                        ('grid_validation', 'passed_full_independent_baliphy_categorical_grid_and_checkpoint_contracts')]:
        record = json.loads(Path(plan[key]).read_text()); assert record['status'] == status
        for path,digest in record['source_hashes'].items(): assert sha(path) == digest
    assert psutil.virtual_memory().available >= 32 * 2**30
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    root = Path(plan['output']); assert not root.exists()
    assert plan['resources']['cpu_per_stage'] == 2 and plan['resources']['memory_gib'] == 32
    assert plan['resources']['swap_gib'] == 0 and plan['resources']['blas_threads'] == 1
    dependencies = [plan['diagnostic_closure_launch'], plan['scalar_closure_launch']]
    for path in dependencies:
        record = json.loads(Path(path).read_text()); record['launch'] = path
        assert fingerprint(record) is None; journal_terminal(record)
    producer = launch('independent-baliphy-category', 'metadata/independent_baliphy_category_wait_plan_20261002.json',
        [sys.executable, 'scripts/prepare_independent_baliphy_category_readback.py', '--plan', str(args.plan)], dependencies, str(args.plan))
    reader = launch('independent-baliphy-category-readback', 'metadata/independent_baliphy_category_readback_wait_plan_20261002.json',
        [sys.executable, 'scripts/readback_independent_baliphy_category_readback.py', '--plan', str(args.plan),
         '--output', str(root / 'readback.json')], [producer], str(args.plan))
    completion = dict(source_plan=str(args.plan), producer_receipt=str(root / 'receipt.json'),
        independent_readback=str(root / 'readback.json'),
        producer_status='complete_full_independent_baliphy_categorical_comparison_pending_readback',
        reader_status='passed_full_independent_baliphy_categorical_serialized_comparison_readback',
        completed_status='complete_verified_full_independent_baliphy_categorical_comparison',
        summary_fields=SUMMARY_FIELDS, launches=[producer, reader],
        pins={p:sha(p) for p in [str(args.plan), producer, reader,
            'scripts/close_full_triad_sequence_stage.py', 'scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'], scope=plan['scope'])
    create(plan['completion_plan'], completion)
    closer = launch('independent-baliphy-category-closure', 'metadata/independent_baliphy_category_closure_wait_plan_20261002.json',
        [sys.executable, 'scripts/close_full_triad_sequence_stage.py', '--plan', plan['completion_plan']],
        [producer, reader], plan['completion_plan'])
    create(plan['launch_inventory'], dict(status='queued_full_independent_baliphy_categorical_numerical_comparison',
        source_plan=str(args.plan), source_plan_sha256=sha(args.plan), launches=[producer, reader, closer],
        native_sampling_restarted=False, gpu=False, new_cost_usd=0, scope=plan['scope']))


if __name__ == '__main__': main()
