#!/usr/bin/env python3
"""Launch all 931 matched-predictor fits after complete input and software proof."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
import sys

import psutil

from ancestral_chain_attempt import sha
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from matched_predictor_branch_inputs import verify
from matched_predictor_branch_fits import load,SUMMARY_FIELDS
from run_matched_predictor_branch_fits import STATUS as PRODUCER_STATUS
from readback_matched_predictor_branch_fits import STATUS as READER_STATUS


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--validation',type=Path,required=True);a=p.parse_args()
    gate=json.loads(a.validation.read_text());assert gate['status']=='passed_full_matched_predictor_native_fit_software_contracts'
    assert (gate['full_real_inputs'],gate['full_real_native_role_configs'],gate['synthetic_native_roles'])==(133,931,7)
    assert gate['independent_dendropy_tree_and_report_readback'] and gate['full_mock_role_serialization_checked']
    assert gate['artificial_failed_roles_retained']==gate['artificial_unresolved_inputs']==2
    assert len(gate['malformed_cases_rejected'])==8;verify(gate['source_hashes'])
    ep=Path('metadata/matched_predictor_branch_fits_software_execution_20261003_v2.json');execution=json.loads(ep.read_text())
    assert execution['status']=='exited_zero_with_receipt' and execution['exit_code']==0 and execution['receipt_sha256']==sha(a.validation)
    verify(execution['source_hashes'])
    native=Path(shutil.which('iqtree3')).resolve();assert str(native) in gate['source_hashes']
    assert gate['source_hashes'][str(native)]==sha(native)
    assert psutil.virtual_memory().available>=12*2**30 and shutil.disk_usage('.').free>=128*2**30
    input_plan='metadata/matched_predictor_branch_inputs_plan_20261003_v1.json'
    input_completion='metadata/matched_predictor_branch_inputs_completed_20261003_v1.json'
    path='metadata/matched_predictor_branch_fits_plan_20261003_v1.json'
    output='results/phylogeny/matched-predictor-branch-fits-20261003-v1';assert not Path(output).exists()
    own=['matched_predictor_branch_fits','run_matched_predictor_branch_fits','readback_matched_predictor_branch_fits',
         'check_matched_predictor_branch_fits','prepare_and_launch_matched_predictor_branch_fits','ancestral_chain_attempt',
         'run_paired_marker_fits','launch_baliphy_reference_sampler_qualification','launch_full_triad_sequence_geometry_followup',
         'close_full_triad_sequence_stage','record_completed_process_handoffs_v2','run_after_verified_dependencies_v2']
    pins={f'scripts/{n}.py':sha(f'scripts/{n}.py') for n in own}
    for q in [str(a.validation),str(ep),input_plan,input_completion,'metadata/3di_substitution_model_receipt.json',
              'metadata/matched_predictor_branch_fits_software_resources_20261003_v1.json',str(native)]:pins[q]=sha(q)
    models={k:str(Path('data/structural_models/garg-hochberg-v3')/('Q.3Di.'+k)) for k in ['AF','LLM']}
    for q in models.values():pins[q]=sha(q)
    scope=('Full matched-predictor control: all931native roles from133exact shared-observation inputs, '
        'covering4970ready marker/tree comparisons among the full8750case/125marker/70view grid. '
        'All3780insufficient cases remain in the separately closed input grid; overlap is21taxa/71markers, '
        'at most16taxa per fit. Exactly oneAA LG+F+G4 andtwo predictors timesAF+G4/AF+F+G4/LLM+G4, '
        'fixed compatible unrooted topology andidentical original position/taxon/unknown masks; allnative '
        'branches initialized at0.1. Predictor counterparts share a deterministic model-specific seed. '
        'Every native attempt has explicit binary/model/input hashes,PID/create/command identity,wall/CPU/AS/file bounds; '
        'one attempt only,failed/invalid outcomes retained,no automaticretry ororphan adoption. Independent '
        'DendroPy tree/report readback plus full source/artifact/two-original-journal closure required. '
        'Point rates are expected state substitutions/site,not physical displacement/year or calibrated '
        'evolutionary acceleration. Native numerical integrity is not accepted model/framework,orthology, '
        'uncertainty,experimental predictor error or aim completion. NoGPU inference/new costs or changes '
        'toexistingancestral/speciestree/covariance jobs. Broader fungal analyses remain required.')
    resources=dict(checked_utc=datetime.now(timezone.utc).isoformat(),cpus=4,workers=4,memory_gib=12,swap_gib=0,
        blas_threads=1,native_threads=1,native_address_space_gib=2,native_cpu_seconds=300,native_wall_seconds=600,
        native_per_file_limit_mib=8,planning_output_gib=96,minimum_free_disk_gib=128,
        native_role_count=931,max_taxa_per_fit=16,max_columns=394,naive_sum_role_wall_ceiling_hours=931*600/3600,
        naive_sum_role_cpu_ceiling_hours=931*300/3600,maximum_address_space_reservations_gib=8,
        runtime_uncalibrated=True,finish_eta=None,gpu=False,new_cost_usd=0,
        available_memory_gib=psutil.virtual_memory().available/2**30,available_disk_gib=shutil.disk_usage('.').free/2**30,
        scope='Native per-process caps and four-worker/12GiB/no-swap cgroup separate from planning; bounds terminate failed roles and do not estimate successful optimization runtime. Conservative96GiB output allowance covers small native trees/reports/logs/receipts; free-spaceguard checked atstart and each completed role. No finish ETA.')
    plan=dict(input_plan=input_plan,input_completion=input_completion,executable=str(native),models=models,pins=pins,
        output=output,completion='metadata/matched_predictor_branch_fits_completed_20261003_v1.json',
        launch_inventory='metadata/matched_predictor_branch_fits_launches_20261003_v1.json',resources=resources,scope=scope)
    # Revalidate the entire closed full input grid before writing a launchable plan.
    qualification=a.validation.parent/'matched_predictor_branch_fits_prelaunch_source_plan_20261003_v1.json'
    create(qualification,plan);source,bindings=load(plan,qualification);assert len(source['configs'])==133
    plan['pins'].update(bindings);create(path,plan)
    producer=launch('matched-predictor-branch-fits',[sys.executable,'scripts/run_matched_predictor_branch_fits.py',
        '--plan',path],[],path,cpus=4,memory=12)
    reader=launch('matched-predictor-branch-fits-readback',[sys.executable,'scripts/readback_matched_predictor_branch_fits.py',
        '--plan',path],[producer],path,cpus=2,memory=12)
    completion_plan='metadata/matched_predictor_branch_fits_completion_plan_20261003_v1.json'
    spec=dict(source_plan=path,producer_receipt=output+'/receipt.json',independent_readback=output+'/readback.json',
        producer_status=PRODUCER_STATUS,reader_status=READER_STATUS,completed_status='complete_verified_full_matched_predictor_native_point_fits',
        summary_fields=SUMMARY_FIELDS,launches=[producer,reader],pins={q:sha(q) for q in [path,producer,reader,
            'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},output=plan['completion'],scope=scope)
    create(completion_plan,spec)
    closer=launch('matched-predictor-branch-fits-closure',[sys.executable,'scripts/close_full_triad_sequence_stage.py',
        '--plan',completion_plan],[producer,reader],completion_plan,cpus=2,memory=12)
    create(plan['launch_inventory'],dict(status='launched_full_matched_predictor_native_point_fits',source_plan=path,
        source_plan_sha256=sha(path),launches=[producer,reader,closer],expected_native_roles=931,gpu=False,new_cost_usd=0))


if __name__=='__main__':main()
