#!/usr/bin/env python3
"""Prepare immutable full-grid plans, then queue timing behind retained closure."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

import psutil

from ancestral_chain_attempt import sha
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from prepare_full_retained_shared_entity_resources import inventory
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from reference_measurement_union_sources import verify


RESOURCE = 'metadata/full_retained_shared_entity_resource_inventory_20261003_v1.json'
FIT = 'metadata/full_retained_shared_entity_fit_draft_plan_20261003_v1.json'
TIMING = 'metadata/full_retained_shared_entity_timing_plan_20261003_v1.json'
DEPENDENCY = 'metadata/reduced_covariance_process_source_v3_closure_launch_20261003.json'
SUMMARY = ['logical_cases', 'model_setting_rows', 'unique_cohorts', 'candidate_rows',
           'eligible_candidates', 'timing_groups', 'source_status_counts', 'timing_status_counts',
           'measured_candidate_coverage', 'unmeasured_review_candidate_coverage',
           'conditional_budget_weighted_seconds']


def gates():
    bindings = {}
    for name, version, session, status in [
        ('fitting', 'v2', 82427,
         'passed_full_retained_shared_entity_source_grid_checkpoint_spectral_contracts_v2'),
        ('timing', 'v1', 63214,
         'passed_full_retained_shared_entity_timing_grid_checkpoint_numerical_probe_contracts_v1'),
    ]:
        prefix = f'metadata/full_retained_{name}_software_'
        gp, ep, tp = [prefix + kind + '_20261003_' + version + '.json'
                      for kind in ['validation', 'execution', 'transport']]
        gate, execution, transport = [json.loads(Path(p).read_text()) for p in [gp, ep, tp]]
        assert gate['status'] == status and gate['scientific_eligibility'] is False
        assert gate['source_and_journal_fixtures_synthetic'] is True
        assert execution['status'] == 'exited_zero_with_receipt' and execution['exit_code'] == 0
        assert execution['receipt_sha256'] == sha(gp) and not execution['timed_out']
        assert transport['tool_session_id'] == session and transport['actual_tool_terminal_exit_code'] == 0
        assert transport['wrapper'] == execution['wrapper']
        assert transport['invocation_id'] == execution['invocation_id']
        assert transport['exact_wrapper_pid_journal_entries'] == 1
        assert transport['original_command_start_records'] == 1
        assert transport['original_completion_resource_records'] == 1
        for mapping in [gate['source_hashes'], execution['source_hashes'],
                        execution['artifacts'], transport['source_hashes']]:
            verify(mapping)
        for p in [gp, ep, tp]:
            bindings[p] = sha(p)
        if name == 'fitting':
            assert (gate['synthetic_declared_cohorts'], gate['synthetic_candidate_rows'],
                    gate['synthetic_setting_fit_links'], gate['positive_serialized_candidate_cases']) == (5, 6000, 12000, 16)
            assert gate['positive_rehashed_identity_cases_rejected'] == 48
            assert len(gate['rehashed_grid_cases_rejected']) == 16
            assert len(gate['rehashed_source_cases_rejected']) == 6
            assert gate['inherited_numerical_envelopes_never_replaced'] is True
        else:
            assert (gate['synthetic_cases'], gate['synthetic_cohorts'], gate['actual_numeric_groups'],
                    gate['actual_numeric_points'], gate['synthetic_full_eligible_selection_census'],
                    gate['synthetic_selection_groups']) == (6000, 5, 8, 16, 1200, 40)
            assert len(gate['rehashed_census_alterations_rejected']) == 6
            assert len(gate['malformed_numeric_probe_exports_rejected']) == 10
            assert gate['source_validation_and_both_producer_guards_and_reader_latent_guard_costs_included'] is True
    return bindings


def prepare():
    bindings = gates()
    for p in [RESOURCE, FIT, TIMING]:
        assert not Path(p).exists(), 'Prepared plans are immutable'
    resources = inventory(); create(RESOURCE, resources)
    original = 'metadata/full_shared_entity_fit_draft_plan_20261002_v2.json'
    fit = json.loads(Path(original).read_text()); verify(fit['pins'])
    reduced_path = 'metadata/full_reduced_covariance_qualification_plan_20261003_v3.json'
    reduced = json.loads(Path(reduced_path).read_text()); verify(reduced['pins'])
    assert reduced['qualification_completion'] == 'metadata/full_uniform_covariance_process_completed_20261003_v2.json'
    fit.update(status='software_qualified_full_retained_draft_pending_real_source_closure_and_timing',
               prepared_utc=datetime.now(timezone.utc).isoformat(),
               qualification_completion=reduced['qualification_completion'],
               retained_plan=reduced_path, retained_completion=reduced['completion'],
               output='results/phylogeny/full-retained-shared-entity-fits-20261003-v1',
               fixture_validation='metadata/full_retained_fitting_software_validation_20261003_v2.json',
               reader_script='scripts/readback_full_retained_shared_entity_fits.py',
               launch_state='not_launched_or_queued', source_original_draft=original)
    assert not Path(fit['output']).exists()
    fit['resources']['resource_inventory'] = RESOURCE
    fit['resources']['runtime_calibrated'] = False
    fit['pins'].update(bindings)
    own = ['full_retained_shared_entity_fit_sources', 'prepare_full_retained_shared_entity_fits',
           'readback_full_retained_shared_entity_fits', 'retained_shared_entity_candidate',
           'readback_retained_shared_entity_candidate', 'full_reduced_covariance_sources_v3',
           'run_full_reduced_covariance_qualification_v3', 'full_retained_shared_entity_timing',
           'prepare_full_retained_shared_entity_timing', 'readback_full_retained_shared_entity_timing',
           'prepare_full_retained_shared_entity_resources',
           'prepare_and_launch_full_retained_shared_entity_timing']
    paths = ['scripts/' + n + '.py' for n in own] + [RESOURCE, original, reduced_path, DEPENDENCY,
        'metadata/process_covariance_readback_plan_20261003_v2.json',
        'metadata/process_covariance_readback_v2_closure_launch_20261003.json',
        'metadata/full_exact_covariance_folds_completed_20261003_v2.json',
        'scripts/launch_baliphy_reference_sampler_qualification.py',
        'scripts/run_after_verified_dependencies_v2.py', '/usr/bin/prlimit', sys.executable]
    fit['pins'].update({p: sha(p) for p in paths}); verify(fit['pins'])
    fit['scope'] = (
        'Complete original4340cohorts/130200designs/260400fit inputs/622080settings: '
        '5208000potential ML/REML candidates and12441600setting-method links, both outcomes/modes '
        'and allfive trees. Same design/response arrays and frozen bounded multistart optimizers. '
        'Closed original seven-kernel process readback, exact nonnegative-cone certificates and '
        'closed retained V3 qualification are mandatory. Every original/retained audit and exact '
        'certificate remains attached. Both q4/q5 bases are selected without removing exception '
        'kernels or altering inherited numerical envelopes; fresh production and independent latent '
        'guards remain additional requirements. The original draft tolerances1e-6 and column_batch32 '
        'are preserved; synthetic numerical checks used gradient3e-6 and smaller column batches. '
        'All unready/constant/source/optimizer/reader review states remain explicit. Restartable '
        'full exports and independent source/spectral/curvature/search/link readback are software '
        'qualified; actual whole-grid timing and enforced resources still gate any fitting launch. '
        'No fitting is launched or queued by this draft. Nonuniform controls, inference/calibration, '
        'accepted phylogenetic/reconciliation/dating framework and all eight biological aims remain required.')
    create(FIT, fit)
    timing = dict(fit_plan=FIT, fit_plan_sha256=sha(FIT), scaled_variance_points=[0.0, 1.0],
        retained_closure_launch=DEPENDENCY,
        output='results/phylogeny/full-retained-shared-entity-input-timing-20261003-v1',
        completion='metadata/full_retained_shared_entity_timing_completed_20261003_v1.json',
        closure_plan='metadata/full_retained_shared_entity_timing_completion_plan_20261003_v1.json',
        launch_inventory='metadata/full_retained_shared_entity_timing_launches_20261003_v1.json',
        resources=resources['timing'], summary_fields=SUMMARY,
        pins={p: sha(p) for p in paths + [FIT, *bindings]},
        scope='Whole original5208000candidate census and4340cohorts, all methods/outcomes/modes/five trees. '
              'One deterministic largest-active-column/highest-condition qualified candidate per eligible '
              'cohort/mode/tree/method/outcome; preserve every exclusion/review and exact identity. Two '
              'variance points0and1 measure source validation, constructors, inherited/fresh production '
              'guards, independent latent qualification and primary/component-spectral agreement. '
              'Checkpoint, complete original source/census/selection/numerical-export/planning readback '
              'and two original journals gate accounting closure. Claimed review groups are remeasured. '
              'Conditional extrapolation is neither a guaranteed bound nor a finish ETA; unmeasured '
              'reviews/loading/export costs remain explicit. No production optimization, altered source, '
              'failed-job restart, GPU work, new charges or biological acceptance.')
    assert not Path(timing['output']).exists(); create(TIMING, timing)
    print(json.dumps(dict(status='prepared_full_retained_fit_and_timing_plans',
                          fit_plan=FIT, timing_plan=TIMING, production_fitting_launched=False)))


def queue():
    gates()
    plan = json.loads(Path(TIMING).read_text()); verify(plan['pins'])
    assert sha(plan['fit_plan']) == plan['fit_plan_sha256']
    fit = json.loads(Path(FIT).read_text()); verify(fit['pins'])
    assert fit['launch_state'] == 'not_launched_or_queued' and not Path(fit['output']).exists()
    for p in [plan['output'], plan['completion'], plan['closure_plan'], plan['launch_inventory']]:
        assert not Path(p).exists()
    assert psutil.virtual_memory().available >= 32 * 2**30
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    original = json.loads(Path(DEPENDENCY).read_text()); original['launch'] = DEPENDENCY
    assert sha(original['plan']) == original['plan_sha256']
    if fingerprint(original) is None:
        journal_terminal(original)
    prefix = ['/usr/bin/prlimit', '--as=' + str(24 * 2**30), '--cpu=604800',
              '--fsize=' + str(2 * 2**30), '--']
    producer = launch('full-retained-shared-entity-timing-v1',
        prefix + [sys.executable, 'scripts/prepare_full_retained_shared_entity_timing.py', '--plan', TIMING],
        [DEPENDENCY], TIMING, cpus=2, memory=32)
    root = Path(plan['output'])
    reader = launch('full-retained-shared-entity-timing-v1-readback',
        prefix + [sys.executable, 'scripts/readback_full_retained_shared_entity_timing.py',
                  '--plan', TIMING, '--output', str(root / 'readback.json')],
        [producer], TIMING, cpus=2, memory=32)
    create(plan['closure_plan'], dict(source_plan=TIMING, producer_receipt=str(root / 'receipt.json'),
        independent_readback=str(root / 'readback.json'),
        producer_status='complete_full_retained_shared_entity_timing_pending_accounting_readback_v1',
        reader_status='passed_full_retained_shared_entity_timing_census_and_planning_accounting_readback_v1',
        completed_status='complete_verified_full_retained_shared_entity_timing_accounting_v1',
        summary_fields=SUMMARY, launches=[producer, reader],
        pins={p: sha(p) for p in [TIMING, producer, reader, 'scripts/close_full_triad_sequence_stage.py',
                                 'scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'], scope=plan['scope']))
    closer = launch('full-retained-shared-entity-timing-v1-closure',
        [sys.executable, 'scripts/close_full_triad_sequence_stage.py', '--plan', plan['closure_plan']],
        [producer, reader], plan['closure_plan'], cpus=2, memory=32)
    create(plan['launch_inventory'], dict(status='queued_full_retained_input_timing_after_retained_source_closure',
        source_plan=TIMING, source_plan_sha256=sha(TIMING), launches=[producer, reader, closer],
        full_potential_candidates=5208000, maximum_eligible_timing_groups=173600,
        production_fitting_launched=False, gpu=False, new_cost_usd=0, all_eight_aims_incomplete=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--action', choices=['prepare', 'queue'], required=True)
    args = parser.parse_args()
    prepare() if args.action == 'prepare' else queue()
