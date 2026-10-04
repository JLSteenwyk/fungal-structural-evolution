#!/usr/bin/env python3
"""Prepare operational fit resources only after unchanged measured inputs close.

This preparation command does not submit a worker. Pending prerequisite
receipts are reported without creating a fitting root or installing resources.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys
import psutil

from ancestral_chain_attempt import sha
from full_weighted_parallel_fit_admission_v1 import closed_timing, production_contract, capacity, FULL
from full_weighted_fit_exports import atomic, PRODUCER, READER, SUMMARY

from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record
from reference_measurement_union_sources import verify


def headers(request):
    fp=Path(request['fit_plan']);tp=Path(request['timing_plan'])
    assert sha(fp)==request['fit_plan_sha256'] and sha(tp)==request['timing_plan_sha256']
    fit=json.loads(fp.read_text());timing=json.loads(tp.read_text())
    assert request['expected']==fit['expected']==FULL
    assert timing['fit_plan']==str(fp) and timing['fit_plan_sha256']==sha(fp)
    assert request['timing_completion']==timing['completion']
    assert fit['launch_state']=='not_launched_or_queued' and not Path(fit['output']).exists()
    production_contract(fit,dict(logical_cases=75188,unique_cohorts=4340,candidate_rows=20832000,model_setting_rows=622080))
    verify(fit['pins']);verify(timing['pins'])
    return fit,timing


def requirements(request):
    fit,timing=headers(request)
    missing=[p for p in [fit['qualification_completion'],timing['completion']] if not Path(p).is_file()]
    return fit,timing,missing


def run(request_path):
    request_path=Path(request_path);request=json.loads(request_path.read_text())
    verify(request['pins']);fit,timing,missing=requirements(request)
    if missing:
        handles=[]
        for path in request['original_closure_launches']:
            r=json.loads(Path(path).read_text());r['launch']=path;assert sha(r['plan'])==r['plan_sha256']
            process=fingerprint(r);handles.append(journal_terminal(r) if process is None else live_record(r,process))
        return dict(status='pending_original_numerical_and_timing_closure_no_fit_prepared',
            checked_utc=datetime.now(timezone.utc).isoformat(),missing_prerequisites=missing,original_handles=handles,
            fit_resources_installed=False,production_fitting_launched_or_queued=False,scientific_eligibility=False,
            request_sha256=sha(request_path),gpu=False,new_cost_usd=0)
    # Successful preparation is separately software gated; waiting on an
    # original job does not install fit resources or imply numerical success.
    gp=Path(request['admission_software_validation']);proof_path=Path(request['admission_software_transport'])
    v=json.loads(gp.read_text());proof=json.loads(proof_path.read_text())
    assert v['status']=='passed_complete_parallel_timing_checkpoint_fit_admission_v1'
    assert v['total_admission_candidate_rows']==72000 and v['full_checkpoint_pairs_checked']==15
    assert len(v['private_rehashed_corruptions_rejected'])==84 and v['synthetic_scope_cannot_admit_production'] is True
    assert proof['original_tool_session_id']==77311 and proof['original_tool_terminal_exit_code']==0
    assert proof['whole_wrapper_initial_and_terminal_payloads_matched']
    assert proof['source_hashes'][str(gp)]==sha(gp)
    verify(v['source_hashes']);verify(proof['source_hashes'])
    export_paths=[Path(request['fitting_software_'+kind]) for kind in ['validation','execution','transport']]
    ev,ee,et=[json.loads(p.read_text()) for p in export_paths]
    assert ev['status']=='passed_complete_guarded_parallel_source_fitting_exports_and_numeric_reader_v3'
    assert ev['total_candidate_rows']==72000 and ev['total_setting_fit_links']==144000
    assert ev['actual_serialized_numeric_replays']==64 and len(ev['private_rehashed_corruptions_rejected'])==23
    assert len(ev['in_memory_mutations_rejected'])==18 and ev['actual_guarded_native_candidate_and_reader_checks']==64
    assert ev['qualified_parent_artifacts_preserved'] and ev['full_export_native_fit_branches_explicitly_mocked']
    assert ee['status']=='exited_zero_with_receipt' and ee['exit_code']==0 and not ee['timed_out']
    assert ee['receipt_sha256']==sha(export_paths[0])
    assert et['original_tool_session_id']==12762 and et['original_tool_terminal_exit_code']==0
    assert et['whole_wrapper_initial_and_terminal_payloads_matched']
    assert et['wrapper']==ee['wrapper'] and et['invocation_id']==ee['invocation_id']
    assert et['manager_start_records']==et['manager_completion_records']==1
    for record in [ev,ee,et]:verify(record['source_hashes']);verify(record.get('artifacts',{}))
    source,bindings,admission=closed_timing(request['fit_plan'],request['timing_plan'],request['timing_completion'])
    production_contract(fit,admission);resources=capacity(admission,fit)
    assert psutil.virtual_memory().available>=resources['minimum_available_ram_gib']*2**30
    assert shutil.disk_usage('.').free>=resources['minimum_free_disk_gib']*2**30
    for path in request['original_closure_launches']:
        r=json.loads(Path(path).read_text());r['launch']=path
        assert fingerprint(r) is None, 'Full original closer must be terminal'
        journal_terminal(r)
    for path in [request['operational_output'],request['resources_output'],request['admission_output']]:assert not Path(path).exists()
    bindings.update(request['pins']);bindings.update({str(p):sha(p) for p in [request_path,gp,proof_path,*export_paths]})
    verify(bindings)
    admission.update(source_hashes=bindings,source_contract=source['fit_contract'],scientific_eligibility=False)
    atomic(Path(request['admission_output']),admission);atomic(Path(request['resources_output']),resources)
    root=Path(fit['output'])
    operation=dict(status='prepared_full_weighted_fit_operations_not_launched',fit_plan=request['fit_plan'],
        fit_plan_sha256=sha(request['fit_plan']),fit_contract=source['fit_contract'],expected=FULL,
        admission=request['admission_output'],resources=request['resources_output'],
        pins={str(p):sha(p) for p in [request_path,request['admission_output'],request['resources_output'],gp,proof_path,
            request['fit_plan'],request['timing_plan'],request['timing_completion'],*export_paths,Path(__file__),Path('scripts/weighted_parallel_fit_execution_runtime_v2.py'),
            Path('scripts/prepare_full_weighted_parallel_source_fits_v2.py'),Path('scripts/readback_full_weighted_parallel_source_fits_v2.py'),
            Path('scripts/weighted_fit_numeric_memory_guard_v1.py')]},
        producer_command=[sys.executable,'scripts/prepare_full_weighted_parallel_source_fits_v2.py','--plan',request['fit_plan']],
        reader_command=[sys.executable,'scripts/readback_full_weighted_parallel_source_fits_v2.py','--plan',request['fit_plan'],'--output',str(root/'readback.json')],
        producer_status=PRODUCER,reader_status=READER,summary_fields=SUMMARY,
        resource_enforcing_controller_implemented=True,parallel_checkpoint_source_adapter_qualified=True,
        cached_and_native_input_memory_guards_required=True,launch_state='not_launched_or_queued',
        fits_computed=0,scientific_eligibility=False,nonuniform_weighting_accepted=False,
        component_variance_attribution_accepted=False,gpu=False,new_cost_usd=0,
        scope='Complete original measured model configuration remains unchanged so source/candidate identities match qualified timing. Operational admission/capacity and command descriptors are separate; no resource installed or native worker submitted. A controller enforcing actual cgroup/BLAS/native caps, original dependencies, all operational pins and complete original fit/readback/journal closure remains required before execution. Scientific calibration and all eight aims remain open.')
    atomic(Path(request['operational_output']),operation)
    return operation


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--request',type=Path,required=True)
    result=run(p.parse_args().request);print(json.dumps(result,indent=2),flush=True)
