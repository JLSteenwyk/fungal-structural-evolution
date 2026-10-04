#!/usr/bin/env python3
"""Gate the complete original source census; no covariance numerical fits."""
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

import psutil

from ancestral_chain_attempt import sha
from full_weighted_covariance_sources_v2 import STATUSES
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from reference_measurement_union_sources import verify
from run_full_weighted_covariance_source_census_v3 import PRODUCER, READER, SUMMARY


def main():
    software = ['metadata/full_weighted_covariance_sources_software_' + k + '_20261004_v4.json'
        for k in ['validation', 'execution', 'transport']]
    v, e, t = [json.loads(Path(p).read_text()) for p in software]
    assert v['status'] == 'passed_complete_declared_four_control_source_census_contracts_v1'
    assert e['status'] == 'exited_zero_with_receipt' and e['exit_code'] == 0 and not e['timed_out']
    assert e['receipt_sha256'] == sha(software[0])
    assert t['actual_tool_session_id'] == 29272 and t['actual_tool_terminal_exit_code'] == 0
    assert t['wrapper'] == e['wrapper'] and t['invocation_id'] == e['invocation_id']
    assert t['exact_wrapper_pid_journal_entries'] == 2 and t['original_start_records'] == t['original_completion_records'] == 1
    assert (v['synthetic_cohorts'], v['synthetic_designs'], v['synthetic_fit_inputs'], v['synthetic_settings'],
        v['prospective_audits'], v['prospective_setting_links']) == (5, 150, 300, 600, 6000, 24000)
    assert len(v['rehashed_output_cases_rejected']) == 8 and len(v['rehashed_source_cases_rejected']) == 10
    assert v['exact_uniform_constant_near_uniform_route_checks'] == 30
    assert v['completed_stage_restart_refusals'] is v['substantive_positive_restoration'] is True
    assert v['source_and_journal_fixtures_synthetic'] is True and v['numerical_audits_computed'] == 0
    for r in [v, e, t]: verify(r['source_hashes']); verify(r.get('artifacts', {}))
    plans = dict(operator='metadata/full_entity_operator_plan_20261002.json',
        design='metadata/full_expanded_model_designs_plan_20261002_v2.json',
        inputs='metadata/full_expanded_model_inputs_plan_20261002.json',
        covariance='metadata/full_expanded_covariance_plan_20261002.json',
        exact='metadata/full_exact_covariance_folds_plan_20261003_v2.json',
        cone='metadata/full_nonuniform_covariance_cone_plan_20261003_v1.json',
        weights='metadata/full_inverse_reuse_weight_plan_20261003_v1.json')
    completions = dict(operator='metadata/full_entity_operator_bank_completed_20261002.json',
        design='metadata/full_expanded_model_designs_v2_completed_20261002.json',
        inputs='metadata/full_expanded_model_inputs_completed_20261002.json',
        covariance='metadata/full_expanded_covariance_completed_20261002.json',
        exact='metadata/full_exact_covariance_folds_completed_20261003_v2.json',
        cone='metadata/full_nonuniform_covariance_cones_completed_20261003_v1.json',
        weights='metadata/full_inverse_reuse_weights_completed_20261003_v1.json')
    parents = {}; archives = []
    for k, p in completions.items():
        c = json.loads(Path(p).read_text()); assert c['status'] == STATUSES[k]
        assert c['scientific_eligibility'] is False and c['exact_process_journals_checked'] == 2
        assert c['logical_cases'] == 75188 and sha(c['full_hash_archive']) == c['full_hash_archive_sha256']
        parents[k] = c; archives.append(c['full_hash_archive'])
    d = parents['design']; w = parents['weights']
    assert (d['unique_cohorts'], d['unique_designs'], d['unique_fit_inputs'], d['model_setting_rows'],
        d['cohort_member_occurrences'], d['largest_cohort']) == (4340, 130200, 260400, 622080, 34110120, 22934)
    assert w['cohorts'] == 4340 and w['case_row_occurrences'] == 34110120
    assert w['residual_diagonal_prepared'] is True and w['raw_reml_basis_qualification_complete'] is False
    resources = dict(checked_utc=datetime.now(timezone.utc).isoformat(), cpus=2, memory_gib=16, swap_gib=0,
        blas_threads=1, address_space_gib=12, cpu_seconds_per_stage=21600, per_file_limit_mib=512,
        minimum_free_disk_gib=128, output_and_scratch_allowance_gib=2,
        logical_cases=75188, original_cohorts=4340, original_designs=130200, original_fit_inputs=260400,
        original_settings=622080, original_cohort_row_occurrences=34110120,
        prospective_numerical_audits=5208000, prospective_setting_audit_links=24883200,
        maximum_cohort_rows=22934, maximum_design_columns=10,
        raw_predictor_bytes_streamed_upper_bound=34110120 * 30 * 10 * 8,
        raw_response_bytes_streamed=34110120 * 30 * 2 * 8,
        source_control_array_bytes=34110120 * 96,
        expanded_species_factor_bytes_upper_bound=75188 * 301 * 5 * 8,
        compressed_cohort_catalog_planning_bytes=4340 * 32768,
        sqlite_scratch_planning_bytes=400 * 2**20,
        runtime_planning_seconds_per_stage=[60, 21600], runtime_calibrated=False,
        wall_time_limit_enforced=False, available_ram_gib=psutil.virtual_memory().available / 2**30,
        free_disk_gib=shutil.disk_usage('.').free / 2**30, new_cost_usd=0, gpu=False,
        scope='Full source/census reconstruction and serialized readback; no covariance Grams or optimizer work. X/y byte estimates are cumulative streamed work, not memory allocation. RAM/AS/CPU/per-file and one-BLAS-thread caps enforced. Runtime range uncalibrated; CPU cap is not ETA or wall limit. Output/scratch allowance is a planning estimate. Full scope retained.')
    assert psutil.virtual_memory().available >= 32 * 2**30 and shutil.disk_usage('.').free >= 128 * 2**30
    rp = 'metadata/full_weighted_covariance_source_census_resources_20261004_v1.json'; create(rp, resources)
    scope = ('Fresh consumed input hashes from seven independently closed parents; all original 4340 cohorts, '
        '75188 cases, 34110120 cohort-row occurrences, 130200 designs, 260400 response inputs and 622080 settings. '
        'Reconstruct raw X and y from all ten original partitions, preserving case/input IDs, active exactly nonzero '
        'columns, original design and response dispositions. Recheck every saved control array/reciprocal recipe '
        'and exact membership identity. Original group assignment rests on the closed independent SQL/Fraction '
        'stage; it is not newly rederived from case-key tables here. Exact named uniform folds only for actual '
        'all-one D, otherwise retain D and target I using the closed positive-diagonal cone, including pair exceptions. '
        'Full original settings checked against unique design/response identities with SQLite; four controls/two modes/'
        'five trees imply 5208000 future numerical audits and 24883200 future links. These are prospective counts, '
        'not computed covariance audits. Serialized census replay plus all consumed source/artifact hashes and '
        'two actual original process journals gate closure. No numerical audit/envelope inherited, weighted fit, '
        'timing calibration, accepted biological effect, posterior acceptance, native restart, GPU or new charge. '
        'Unconsumed broader parent archive artifacts are not freshly replayed. All eight aims remain incomplete.')
    modules = ['full_weighted_covariance_sources_v2', 'run_full_weighted_covariance_source_census_v3',
        'check_full_weighted_covariance_sources_v4', 'prepare_and_launch_full_weighted_covariance_source_census',
        'full_covariance_qualification_sources', 'full_exact_covariance_sources', 'full_expanded_model_design_sources',
        'full_expanded_model_input_sources', 'inverse_reuse_weight_controls', 'nonuniform_covariance_cone',
        'reduced_covariance_basis', 'covariance_exact_folds_v2', 'reference_measurement_union_sources',
        'run_full_inverse_reuse_weights', 'ancestral_chain_attempt', 'launch_baliphy_reference_sampler_qualification',
        'launch_full_triad_sequence_geometry_followup', 'run_after_verified_dependencies_v2',
        'close_full_triad_sequence_stage', 'record_completed_process_handoffs_v2']
    pins = ['scripts/' + m + '.py' for m in modules] + [rp, *software, *plans.values(), *completions.values(),
        *archives, sys.executable, '/usr/bin/prlimit']
    output = 'results/phylogeny/full-four-control-covariance-source-census-20261004-v1'; assert not Path(output).exists()
    pp = 'metadata/full_weighted_covariance_source_census_plan_20261004_v1.json'
    completion = 'metadata/full_weighted_covariance_source_census_completed_20261004_v1.json'
    plan = {k + '_plan': p for k, p in plans.items()}; plan.update({k + '_completion': p for k, p in completions.items()})
    plan.update(trees=d['trees'], pins={p: sha(p) for p in pins}, output=output, resources=resources, scope=scope,
        expected=dict(logical_cases=75188, cohorts=4340, designs=130200, fit_inputs=260400, settings=622080,
            case_row_occurrences=34110120, numerical_audit_rows=5208000, setting_audit_links=24883200), completion=completion)
    create(pp, plan)
    prefix = ['/usr/bin/prlimit', '--as=' + str(12 * 2**30), '--cpu=21600', '--fsize=' + str(512 * 2**20), '--']
    # The existing immutable launcher uses its original date in handle filenames;
    # actual launch UTC, source plan date and exact identities are recorded.
    producer = launch('full-weighted-covariance-source-census-v1', prefix + [sys.executable,
        'scripts/run_full_weighted_covariance_source_census_v3.py', '--plan', pp], [], pp, cpus=2, memory=16)
    reader = launch('full-weighted-covariance-source-census-v1-readback', prefix + [sys.executable,
        'scripts/run_full_weighted_covariance_source_census_v3.py', '--plan', pp, '--reader'], [producer], pp, cpus=2, memory=16)
    cp = 'metadata/full_weighted_covariance_source_census_completion_plan_20261004_v1.json'
    create(cp, dict(source_plan=pp, producer_receipt=output + '/receipt.json', independent_readback=output + '/readback.json',
        producer_status=PRODUCER, reader_status=READER, completed_status='complete_verified_full_four_control_covariance_source_census_v1',
        summary_fields=SUMMARY, launches=[producer, reader], output=completion, scope=scope,
        pins={p: sha(p) for p in [pp, producer, reader, 'scripts/close_full_triad_sequence_stage.py', 'scripts/record_completed_process_handoffs_v2.py']}))
    closer = launch('full-weighted-covariance-source-census-v1-closure', [sys.executable,
        'scripts/close_full_triad_sequence_stage.py', '--plan', cp], [producer, reader], cp, cpus=2, memory=16)
    create('metadata/full_weighted_covariance_source_census_launches_20261004_v1.json',
        dict(status='queued_full_original_four_control_source_census_and_serialized_readback', source_plan=pp,
            source_plan_sha256=sha(pp), launches=[producer, reader, closer], expected=plan['expected'],
            numerical_audits_computed=0, weighted_fitting_launched=False, gpu=False, new_cost_usd=0, all_eight_aims_incomplete=True))


if __name__ == '__main__': main()
