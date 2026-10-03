#!/usr/bin/env python3
"""Queue the complete retained-kernel grid behind original arithmetic closure."""
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

import psutil

from ancestral_chain_attempt import sha
from full_reduced_covariance_sources import SUMMARY_FIELDS
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from reference_measurement_union_sources import verify
from run_full_reduced_covariance_qualification import PRODUCER_STATUS, READER_STATUS


def main():
    gate_path = Path('metadata/reduced_covariance_qualification_software_validation_20261003_v1.json')
    execution_path = Path('metadata/reduced_covariance_qualification_software_execution_20261003_v1.json')
    transport_path = Path('metadata/reduced_covariance_qualification_software_transport_20261003_v1.json')
    gate, execution, transport = [json.loads(p.read_text()) for p in [gate_path, execution_path, transport_path]]
    assert gate['status'] == 'passed_full_exact_retained_covariance_qualification_software_contracts'
    assert (gate['synthetic_pipeline_audits'], gate['synthetic_pipeline_setting_links']) == (600, 1200)
    assert (gate['full_real_certificates_checked'], gate['full_real_cohorts'], gate['full_real_logical_cases']) == (8680, 4340, 75188)
    assert len(gate['dense_and_review_scenarios']) == 7 and len(gate['malformed_cases_rejected']) == 16
    assert gate['dense_nonnegative_cone_roundtrips'] == 80
    for k in ['all_five_trees_and_both_modes', 'inherited_numerical_envelopes_unchanged',
              'nonready_setting_states_retained', 'completed_restart_refused', 'source_and_journal_fixtures_synthetic']:
        assert gate[k] is True
    assert execution['status'] == 'exited_zero_with_receipt' and execution['exit_code'] == 0
    assert execution['receipt_sha256'] == sha(gate_path)
    assert transport['status'] == 'verified_original_retained_covariance_software_wait_exited_zero'
    assert transport['actual_tool_terminal_exit_code'] == 0 and transport['invocation_id'] == execution['invocation_id']
    for r in [gate, execution, transport]:
        verify(r['source_hashes'])
        if 'artifacts' in r: verify(r['artifacts'])
    exact_completion = 'metadata/full_exact_covariance_folds_completed_20261003_v2.json'
    exact = json.loads(Path(exact_completion).read_text())
    assert exact['status'] == 'complete_verified_full_exact_uniform_covariance_folds_v2'
    assert exact['cohorts'] == 4340 and exact['certificates'] == 8680
    assert exact['retained_basis_counts'] == {'4': 7232, '5': 1448}
    assert sha(exact['full_hash_archive']) == exact['full_hash_archive_sha256']
    dependencies = ['metadata/full_uniform_covariance_qualification_closure_launch_20261002.json',
                    'metadata/full_exact_covariance_general_folds_closure_launch_20261003.json']
    observations = []
    for p in dependencies:
        row = json.loads(Path(p).read_text()); row['launch'] = p
        assert sha(row['plan']) == row['plan_sha256']
        process = fingerprint(row)
        if process is None: observations.append(journal_terminal(row))
        else: observations.append(dict(launch=p, pid=process.pid, created=process.create_time(),
            cmdline=process.cmdline(), status=process.status(), scope='Exact original dependency confirmed live; not complete.'))
    assert psutil.virtual_memory().available >= 16 * 2**30 and shutil.disk_usage('.').free >= 64 * 2**30
    output = 'results/phylogeny/full-exact-retained-covariance-qualification-20261003-v1'
    assert not Path(output).exists()
    expected = dict(logical_cases=75188, model_setting_rows=622080, unique_cohorts=4340,
        unique_designs=130200, audit_rows=1302000, setting_audit_links=6220800)
    scope = ('Complete numerical qualification of every exact retained uniform covariance kernel across '
        'all75188 original cases,4340cohorts,130200designs,both loading modes and five working species trees; '
        '1302000audits and6220800original setting links,including every outcome/fit-input identity and review/failure state. '
        'Requires original full seven-kernel arithmetic/latent/SQL/source/journal closure and full exact V2 named '
        'integer operator/cone closure before consumption. Retained four/five kernels are actual unchanged original '
        'kernels; raw/REML Grams and error envelopes are exact principal submatrices. No narrowed numerical bounds, '
        'rank-based deletion,clipping,case filtering or promotion of missing/failed source audits. Independent '
        'entry-wise selection and SciPy gesvd check every output; SQLite identity constraints and lockstep source '
        'links retain all settings. Relevant closed artifacts and full inherited archive hashes are freshly verified; '
        'broader inherited multi-million-binding source files are not rehashed here. Separate component variance '
        'attribution,nonuniform weighting,fitting,optimization/calibration,phylogenetic/model acceptance and all '
        'eight biological aims remain open. No biological pilot,GPU,new charges,old job/source changes or retries.')
    resources = dict(checked_utc=datetime.now(timezone.utc).isoformat(), cpus=2, memory_gib=16, swap_gib=0,
        workers=1, blas_threads=1, address_space_gib=12, cpu_seconds_per_stage=21600, per_file_limit_mib=1024,
        planning_output_gib=8, minimum_free_disk_gib=64, software_child_cpu_seconds=execution['child_cpu_seconds'],
        software_child_peak_rss_bytes=execution['child_peak_rss_bytes'], expected_audits=1302000,
        expected_setting_links=6220800, existing_compressed_audit_bytes=Path(
            'results/phylogeny/full-uniform-covariance-qualification-20261002-v1/design_covariance_audits.jsonl.gz').stat().st_size,
        existing_compressed_link_bytes=Path(
            'results/phylogeny/full-uniform-covariance-qualification-20261002-v1/setting_audit_links.tsv.gz').stat().st_size,
        available_memory_gib=psutil.virtual_memory().available / 2**30, available_disk_gib=shutil.disk_usage('.').free / 2**30,
        finish_eta=None, new_cost_usd=0, gpu=False,
        scope='Planning8GiB output including temporary SQLite index; not an enforced total-byte limit. Original compressed input sizes and software resources are measurements,not full-run estimates. No calibrated ETA. One worker,2CPU/16GiB/no-swap service,12GiB AS/21600CPU-second/1GiB per-file native limits. Dependency waiting precedes native process caps; no automatic retry or fit launch.')
    own = ['prepare_and_launch_full_reduced_covariance_qualification', 'full_reduced_covariance_sources',
        'run_full_reduced_covariance_qualification', 'reduced_covariance_basis',
        'check_full_reduced_covariance_qualification', 'run_reduced_covariance_software_stage_v1',
        'covariance_exact_folds_v2', 'full_exact_covariance_sources', 'covariance_basis_audit',
        'full_covariance_qualification_sources', 'full_expanded_model_design_sources', 'full_expanded_model_input_sources',
        'full_entity_operator_sources', 'background_measurement_union_sources', 'reference_measurement_union_sources',
        'ancestral_chain_attempt', 'run_ortholog_pair_guide_comparison', 'launch_baliphy_reference_sampler_qualification',
        'launch_full_triad_sequence_geometry_followup', 'run_after_verified_dependencies_v2',
        'record_project_runtime_checkpoint_v4', 'close_full_triad_sequence_stage', 'record_completed_process_handoffs_v2']
    pins = {'scripts/' + name + '.py': sha('scripts/' + name + '.py') for name in own}
    for p in [gate_path, execution_path, transport_path, Path('/usr/bin/prlimit'), Path(sys.executable),
        Path(exact_completion), Path('metadata/reduced_covariance_qualification_software_resources_20261003_v1.json'),
        Path('metadata/full_uniform_covariance_qualification_plan_20261002.json'),
        Path('metadata/full_exact_covariance_folds_plan_20261003_v2.json'), *map(Path, dependencies)]: pins[str(p)] = sha(p)
    pins.update(transport['source_hashes']); pins.update(execution['artifacts'])
    path = 'metadata/full_reduced_covariance_qualification_plan_20261003_v1.json'
    plan = dict(qualification_plan='metadata/full_uniform_covariance_qualification_plan_20261002.json',
        qualification_completion='metadata/full_uniform_covariance_qualification_completed_20261002.json',
        exact_plan='metadata/full_exact_covariance_folds_plan_20261003_v2.json', exact_completion=exact_completion,
        dependencies=dependencies, expected=expected, trees=['mafft_guide','pmsf_mafft_profile',
            'pmsf_profile_mafft','pmsf_profile_profile','profile_guide'], output=output, pins=pins, resources=resources,
        completion='metadata/full_reduced_covariance_qualification_completed_20261003_v1.json',
        launch_inventory='metadata/full_reduced_covariance_qualification_launches_20261003_v1.json', scope=scope)
    create(path, plan)
    prefix = ['/usr/bin/prlimit', '--as=' + str(12 * 2**30), '--cpu=21600', '--fsize=' + str(1024 * 2**20), '--']
    producer = launch('full-reduced-covariance-qualification', prefix + [sys.executable,
        'scripts/run_full_reduced_covariance_qualification.py','--plan',path],dependencies,path,cpus=2,memory=16)
    reader = launch('full-reduced-covariance-qualification-readback', prefix + [sys.executable,
        'scripts/run_full_reduced_covariance_qualification.py','--plan',path,'--reader'],[producer],path,cpus=2,memory=16)
    cp = 'metadata/full_reduced_covariance_qualification_completion_plan_20261003_v1.json'
    create(cp, dict(source_plan=path, producer_receipt=output+'/receipt.json', independent_readback=output+'/readback.json',
        producer_status=PRODUCER_STATUS, reader_status=READER_STATUS,
        completed_status='complete_verified_full_exact_retained_uniform_covariance_qualification',
        summary_fields=SUMMARY_FIELDS, launches=[producer,reader], pins={p:sha(p) for p in [path,producer,reader,
            'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'],scope=scope))
    closer = launch('full-reduced-covariance-qualification-closure',[sys.executable,
        'scripts/close_full_triad_sequence_stage.py','--plan',cp],[producer,reader],cp,cpus=2,memory=16)
    create(plan['launch_inventory'],dict(status='queued_complete_exact_retained_covariance_qualification',
        source_plan=path,source_plan_sha256=sha(path),launches=[producer,reader,closer],
        expected=expected,dependency_observations=observations,production_fitting_launched=False,
        raw_reml_basis_qualification_complete=False,gpu=False,new_cost_usd=0,all_eight_aims_incomplete=True))


if __name__ == '__main__': main()
