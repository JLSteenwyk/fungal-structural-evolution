#!/usr/bin/env python3
"""Queue full reread of the original startup outputs; never rerun native jobs."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

from ancestral_chain_attempt import sha
from launch_baliphy_reference_preflight import launch
from launch_full_triad_sequence_geometry_followup import create
from run_baliphy_reference_preflight import SUMMARY_FIELDS, verify


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--native-validation', type=Path, required=True)
    p.add_argument('--workflow-validation', type=Path, required=True); a = p.parse_args()
    native_path = Path('metadata/baliphy_reference_preflight_plan_20261003.json'); native = json.loads(native_path.read_text())
    verify(native)
    for path, status in [(a.native_validation, 'passed_actual_native_reference_startup_footer_contracts'),
                         (a.workflow_validation, 'passed_full_reference_footer_serialization_software_contracts')]:
        r = json.loads(path.read_text()); assert r['status'] == status
        for q, h in r['source_hashes'].items(): assert sha(q) == h, q
    validation = json.loads(a.native_validation.read_text())
    assert validation['actual_full_grid_startups_checked'] >= 16 and validation['actual_timing_footer_records_checked'] > 0
    assert len(validation['malformed_stdout_rejected']) == 7
    workflow = json.loads(a.workflow_validation.read_text())
    assert workflow['full_chains'] == 1620 and workflow['full_quartets'] == 405
    assert len(workflow['altered_exports_and_completed_restarts_rejected']) == 10
    scope = ('All1620startup identities reread from unchanged native attempts after complete v1 provenance closure. '
             'Accept one initial JSON model plus whitespace or exact native five-line Work timing footer; '
             'unexpected output fails. All original dispositions, failures, aliases and seeds retained. '
             'Native-vs-reference homology, free ancestral representation, original model source and degree-aware '
             'density rechecked. No native reruns, posterior sampling, GPU, memory repair or scientific acceptance. '
             'Reader shares the new decoder; full source/artifact/two original journals required for closure.')
    pins = dict(native['pins']); pins[str(native_path)] = sha(native_path)
    paths = [Path(__file__), Path('scripts/baliphy_reference_startup_readback_v2.py'),
             Path('scripts/check_baliphy_reference_startup_footer.py'), Path('scripts/check_baliphy_reference_footer_workflow.py'),
             Path('scripts/launch_baliphy_reference_preflight.py'), Path('scripts/launch_full_triad_sequence_geometry_followup.py'),
             a.native_validation, a.workflow_validation]
    pins.update({str(x): sha(x) for x in paths})
    plan_path = Path('metadata/baliphy_reference_startup_footer_plan_20261003.json')
    output = 'results/ancestral/full-baliphy-reference-startup-footer-replay-20261003-v2'
    source = dict(native_plan=str(native_path), output=output, pins=pins,
        resources=dict(checked_utc=datetime.now(timezone.utc).isoformat(), cpus=2, memory_gib=24,
            swap_gib=0, blas_threads=1, output_allowance_gib=1, minimum_free_disk_gib=64,
            planning_active_hours_per_stage=[.02, 2], runtime_uncalibrated=True,
            finish_eta=None, native_outputs_parsed=1620, new_native_jobs=0, gpu=False, new_cost_usd=0,
            scope='Existing startup JSON/model/receipt/input reread and full archive hashes. Runtime/output/workspace are planning allowances, not calibrated guarantees; no native inference or extra sampling.'),
        completion='metadata/baliphy_reference_startup_footer_completed_20261003.json', scope=scope)
    create(plan_path, source)
    dependency = 'metadata/baliphy_reference_preflight_closure_launch_20261003.json'
    producer = launch('baliphy-reference-startup-footer', [sys.executable, 'scripts/baliphy_reference_startup_readback_v2.py',
        '--plan', str(plan_path)], [dependency], str(plan_path))
    reader = launch('baliphy-reference-startup-footer-readback', [sys.executable, 'scripts/baliphy_reference_startup_readback_v2.py',
        '--plan', str(plan_path), '--reader'], [producer], str(plan_path))
    completion_path = Path('metadata/baliphy_reference_startup_footer_completion_plan_20261003.json')
    completion = dict(source_plan=str(plan_path), producer_receipt=output + '/receipt.json',
        independent_readback=output + '/readback.json',
        producer_status='complete_full_reference_startup_footer_replay_pending_readback',
        reader_status='passed_full_reference_startup_footer_serialized_readback',
        completed_status='complete_verified_full_reference_startup_footer_replay',
        summary_fields=SUMMARY_FIELDS + ['native_timing_footers_checked', 'v1_parser_dispositions_reclassified'],
        launches=[producer, reader], pins={x: sha(x) for x in [str(plan_path), producer, reader,
        'scripts/close_full_triad_sequence_stage.py', 'scripts/record_completed_process_handoffs_v2.py']},
        output=source['completion'], scope=scope)
    create(completion_path, completion)
    closer = launch('baliphy-reference-startup-footer-closure', [sys.executable, 'scripts/close_full_triad_sequence_stage.py',
        '--plan', str(completion_path)], [producer, reader], str(completion_path))
    create('metadata/baliphy_reference_startup_footer_launches_20261003.json',
        dict(status='queued_complete_original_startup_footer_replay', source_plan=str(plan_path),
             source_plan_sha256=sha(plan_path), launches=[producer, reader, closer],
             new_native_jobs=0, posterior_sampling_launched=False, gpu=False, new_cost_usd=0))


if __name__ == '__main__': main()
