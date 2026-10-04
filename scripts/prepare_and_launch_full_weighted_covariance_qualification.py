#!/usr/bin/env python3
"""Gate and queue every four-control numerical audit behind full source closure."""
import csv
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import shutil
import sys
from collections import Counter

import psutil

from ancestral_chain_attempt import sha
from full_weighted_covariance_qualification import PRODUCER, READER, SUMMARY
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from reference_measurement_union_sources import verify


def software(stem, suffix, expected_session):
    paths = ['metadata/' + stem + '_software_' + k + '_' + suffix + '.json' for k in ['validation','execution','transport']]
    v,e,t = [json.loads(Path(p).read_text()) for p in paths]
    assert e['status'] == 'exited_zero_with_receipt' and e['exit_code'] == 0 and not e['timed_out']
    assert sha(paths[0]) == e['receipt_sha256'] and t['actual_tool_terminal_exit_code'] == 0
    assert t['actual_tool_session_id'] == expected_session and t['wrapper'] == e['wrapper'] and t['invocation_id'] == e['invocation_id']
    assert t['exact_wrapper_pid_journal_entries'] == 2 and t['original_start_records'] == t['original_completion_records'] == 1
    for r in [v,e,t]: verify(r['source_hashes']); verify(r.get('artifacts',{}))
    return v,paths


def main():
    primary,ppaths = software('positive_diagonal_basis_context','20261003_v1',34539)
    latent,lpaths = software('independent_positive_diagonal_basis','20261004_v1',84241)
    v,npaths = software('full_weighted_covariance_qualification','20261004_v2',33138)
    assert primary['status'] == 'passed_reused_fresh_positive_diagonal_raw_reml_gram_contracts_v1'
    assert latent['status'] == 'passed_independent_cached_positive_diagonal_basis_contracts_v1'
    assert (latent['general_diagonal_cases'],latent['named_uniform_basis_cases'],latent['row_permutations'],
        latent['full_rank_design_transforms'],latent['reordered_calls'],latent['invalid_diagonal_design_cases_rejected']) == (24,4,24,24,24,40)
    assert v['status'] == 'passed_complete_declared_four_control_covariance_numerical_contracts_v2'
    assert (v['total_declared_audit_records'],v['total_declared_setting_links']) == (18000,72000)
    assert len(v['additional_complete_qualified_grids']) == 2
    assert v['qualified_basis_dimensions'] == [4,5,6]
    assert set(v['qualified_nonuniform_policies']) == {'background_node','background_pair','family_component'}
    assert all(r['qualified_records'] == 1200 and r['qualified_nonuniform_records'] == 900 for r in v['additional_complete_qualified_grids'])
    assert len(v['rehashed_audit_corruptions_rejected']) == 14 and len(v['rehashed_link_corruptions_rejected']) == 5
    assert v['original_numeric_failure_captures'] == 13 and v['completed_stage_restart_refusals'] is v['substantive_positive_restoration'] is True
    assert v['source_and_journal_fixtures_synthetic'] is True and v['scientific_eligibility'] is False
    assert v['unchanged_comparison_rtol'] == 3e-9 and v['unchanged_comparison_atol'] == 2e-8
    sp = Path('metadata/full_weighted_covariance_source_census_plan_20261004_v1.json'); source = json.loads(sp.read_text())
    verify(source['pins'])
    inventory_path = Path('metadata/full_weighted_covariance_source_census_launches_20261004_v1.json')
    inventory = json.loads(inventory_path.read_text()); assert inventory['source_plan_sha256'] == sha(sp)
    source_closer = inventory['launches'][2]; handle = json.loads(Path(source_closer).read_text()); handle['launch'] = source_closer
    if fingerprint(handle) is None: journal_terminal(handle)
    expected = source['expected']; assert expected == dict(logical_cases=75188,cohorts=4340,designs=130200,
        fit_inputs=260400,settings=622080,case_row_occurrences=34110120,numerical_audit_rows=5208000,setting_audit_links=24883200)
    # Complete setting census before any numerical launch; not a fitted-data pilot.
    design_root = Path(json.loads(Path(source['design_plan']).read_text())['output'])
    settings = Counter(); longest = 0
    with gzip.open(design_root / 'model_settings.tsv.gz','rt') as f:
        reader = csv.DictReader(f,delimiter='\t')
        for r in reader:
            settings[r['cohort_id']] += 1; longest = max(longest,len(('\t'.join(r.values())+'\n').encode()))
    assert len(settings) == 4340 and sum(settings.values()) == 622080 and max(settings.values()) == 960 and longest == 348
    available = psutil.virtual_memory().available; disk = shutil.disk_usage('.').free
    assert available >= 64 * 2**30 and disk >= 256 * 2**30
    resources = dict(checked_utc=datetime.now(timezone.utc).isoformat(),cpus=2,memory_gib=32,swap_gib=0,blas_threads=1,
        address_space_gib=24,cpu_seconds_per_stage=604800,per_file_limit_mib=512,
        minimum_free_disk_gib=256,output_and_scratch_allowance_gib=160,
        original_cases=75188,original_cohorts=4340,original_designs=130200,original_settings=622080,
        numerical_audit_records=5208000,setting_audit_links=24883200,
        audits_per_cohort=1200,maximum_original_settings_per_cohort=960,maximum_links_per_cohort=38400,
        measured_longest_original_setting_bytes=348,audit_record_budget_bytes=16384,link_record_budget_bytes=2048,
        uncompressed_audit_budget_bytes=5208000 * 16384,uncompressed_link_budget_bytes=24883200 * 2048,
        maximum_uncompressed_cohort_audit_bytes=1200 * 16384,
        maximum_uncompressed_cohort_link_bytes=960 * 40 * 2048,
        maximum_numeric_failure_input_bytes=22934 * (301 + 10 + 2) * 8,
        output_segment_files=8680,one_cohort_and_mode_tree_context_at_a_time=True,
        runtime_planning_core_hours_per_stage=[6,168],runtime_calibrated=False,wall_time_limit_enforced=False,
        native_cpu_cap_is_eta=False,available_ram_gib=available / 2**30,free_disk_gib=disk / 2**30,
        gpu=False,new_cost_usd=0,
        scope='Complete five-million-record numerical producer and separate latent reader, no optimizer/native structure inference. All original setting counts measured; individual cohort outputs and failure inputs fit512MiB per-file planning bounds. Compressed files should be smaller than conservative raw budgets, but160GiB is a planning allowance, not a directory quota. Native CPU/AS/cgroup/per-file/BLAS caps enforced. Runtime planning range is uncalibrated, not an ETA; no enforced wall limit. Exact full source closure is a hard execution gate.')
    rp = 'metadata/full_weighted_covariance_qualification_resources_20261004_v1.json'; create(rp,resources)
    scope = ('All original4340cohorts/75188cases/130200designs/622080settings, four distinct control labels, '
        'both loading modes/five trees;5208000audit records and24883200original setting links. '
        'Execution gated by exact original full source census closure. Fresh component/image-product primary '
        'and separate long-double/pivoted-QR/gesvd/latent-overlap reader; rebuild every D correction and numerical '
        'envelope. Named source certificates alone select uniform/nonuniform bases; keep targetI for nonuniformD '
        'and preserve genuine physical-pair exceptions, zero norms, dependencies, original nonfit states and '
        'rank boundaries. No tolerance widening, clipping, saved uniform audit/envelope transfer or numerical '
        'basis deletion. Cohort segments retain every audit/setting ordinal and original row fields. Any rejected '
        'numeric case preserves exact inputs/primary/reference/error evidence. Complete every raw/projected '
        'Gram/envelope/diagnostic check and all original link rows, all consumed source/output hashes and two '
        'actual original process journals before closure. No covariance fit, component attribution, weighting '
        'precision calibration, accepted biological effect, posterior qualification or repair of historical '
        'uniform discrepancy implied. No native restart, GPU, new charge or unrelated settings changed. All8aims remain open.')
    modules = ['full_weighted_covariance_qualification','check_full_weighted_covariance_qualification_v2',
        'weighted_covariance_qualified_source_fixture','prepare_and_launch_full_weighted_covariance_qualification',
        'independent_positive_diagonal_basis_context','independent_positive_diagonal_kernel_products',
        'positive_diagonal_basis_context','covariance_basis_context','covariance_basis_audit',
        'readback_full_covariance_qualification','covariance_basis_independent','reduced_covariance_basis',
        'full_weighted_covariance_sources_v2','full_exact_covariance_sources','full_covariance_qualification_sources',
        'full_expanded_model_design_sources','full_expanded_model_input_sources','full_entity_operator_sources',
        'nonuniform_covariance_cone','covariance_exact_folds_v2','inverse_reuse_weight_controls',
        'run_full_inverse_reuse_weights','reference_measurement_union_sources','ancestral_chain_attempt',
        'launch_baliphy_reference_sampler_qualification','launch_full_triad_sequence_geometry_followup',
        'run_after_verified_dependencies_v2','close_full_triad_sequence_stage','record_completed_process_handoffs_v2']
    proof_paths = [*ppaths,*lpaths,*npaths,str(sp),str(inventory_path),source_closer,rp,sys.executable,'/usr/bin/prlimit']
    pins = {p:sha(p) for p in ['scripts/' + m + '.py' for m in modules] + proof_paths}
    output = 'results/phylogeny/full-four-control-covariance-numerical-qualification-20261004-v1'
    assert not Path(output).exists()
    pp = 'metadata/full_weighted_covariance_qualification_plan_20261004_v1.json'
    completion = 'metadata/full_weighted_covariance_qualification_completed_20261004_v1.json'
    plan = {k:v for k,v in source.items() if k.endswith('_plan') or k.endswith('_completion')}
    plan.update(source_census_plan=str(sp),source_census_completion=source['completion'],
        trees=source['trees'],expected=expected,output=output,resources=resources,pins=pins,scope=scope,completion=completion)
    create(pp,plan)
    prefix = ['/usr/bin/prlimit','--as=' + str(24 * 2**30),'--cpu=604800','--fsize=' + str(512 * 2**20),'--']
    producer = launch('full-weighted-covariance-numerical-v1',prefix + [sys.executable,
        'scripts/full_weighted_covariance_qualification.py','--plan',pp],[source_closer],pp,cpus=2,memory=32)
    reader = launch('full-weighted-covariance-numerical-v1-readback',prefix + [sys.executable,
        'scripts/full_weighted_covariance_qualification.py','--plan',pp,'--reader'],[producer],pp,cpus=2,memory=32)
    cp = 'metadata/full_weighted_covariance_qualification_completion_plan_20261004_v1.json'
    create(cp,dict(source_plan=pp,producer_receipt=output + '/receipt.json',independent_readback=output + '/readback.json',
        producer_status=PRODUCER,reader_status=READER,completed_status='complete_verified_full_four_control_covariance_numerical_qualification_v1',
        summary_fields=SUMMARY,launches=[producer,reader],output=completion,scope=scope,
        pins={p:sha(p) for p in [pp,producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']}))
    closer = launch('full-weighted-covariance-numerical-v1-closure',[sys.executable,
        'scripts/close_full_triad_sequence_stage.py','--plan',cp],[producer,reader],cp,cpus=2,memory=32)
    create('metadata/full_weighted_covariance_qualification_launches_20261004_v1.json',
        dict(status='queued_full_original_four_control_numerical_qualification_and_latent_readback',
            source_plan=pp,source_plan_sha256=sha(pp),launches=[producer,reader,closer],expected=expected,
            prerequisite_source_closure=source['completion'],weighted_working_model_fits_launched=False,
            gpu=False,new_cost_usd=0,all_eight_aims_incomplete=True))


if __name__ == '__main__': main()
