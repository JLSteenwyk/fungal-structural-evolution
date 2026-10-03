#!/usr/bin/env python3
"""Launch the complete matched-predictor input grid after bounded software proof."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

import psutil

from audit_selected_taxon_identity_snapshot_v2 import sha
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from matched_predictor_branch_inputs import load, verify, SCHEMA, SUMMARY_FIELDS
from prepare_matched_predictor_branch_inputs import STATUS as PRODUCER_STATUS
from readback_matched_predictor_branch_inputs import STATUS as READER_STATUS


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--validation',type=Path,required=True);a=p.parse_args()
    gate=json.loads(a.validation.read_text());assert gate['status']=='passed_full_matched_predictor_branch_input_software_contracts'
    assert gate['synthetic_full_serialized_comparison_cases']==8750 and gate['independent_raw_tree_pruning_checked']
    assert len(gate['malformed_cases_rejected'])==17
    assert (gate['full_real_tree_views'],gate['full_real_marker_slots'],gate['full_real_taxon_entries'],
            gate['exact_complete_protein_context_cells'])==(70,125,526,673)
    verify(gate['source_hashes'])
    software_execution=Path('metadata/matched_predictor_branch_inputs_software_execution_20261003_v1.json')
    execution=json.loads(software_execution.read_text())
    assert execution['status']=='exited_zero_with_receipt' and execution['exit_code']==0
    assert execution['receipt_sha256']==sha(a.validation);verify(execution['source_hashes'])
    assert psutil.virtual_memory().available>=8*2**30 and shutil.disk_usage('.').free>=16*2**30
    source, bindings=load()
    own=['matched_predictor_branch_inputs','prepare_matched_predictor_branch_inputs','readback_matched_predictor_branch_inputs',
         'check_matched_predictor_branch_inputs','prepare_and_launch_matched_predictor_branch_inputs',
         'launch_baliphy_reference_sampler_qualification','launch_full_triad_sequence_geometry_followup',
         'run_after_verified_dependencies_v2','close_full_triad_sequence_stage','record_completed_process_handoffs_v2']
    paths=[f'scripts/{n}.py' for n in own]+[str(a.validation),str(software_execution)]
    paths+=['metadata/matched_predictor_branch_inputs_software_resources_20261003_v1.json',
            'metadata/structural_marker_tree_projection_completed_20261003_v1.json']
    bindings.update({path:sha(path) for path in paths});verify(bindings)
    plan_path='metadata/matched_predictor_branch_inputs_plan_20261003_v1.json'
    output='results/phylogeny/matched-predictor-branch-inputs-20261003-v1'
    inventory='metadata/matched_predictor_branch_inputs_launches_20261003_v1.json'
    completion='metadata/matched_predictor_branch_inputs_completed_20261003_v1.json'
    assert not Path(output).exists()
    scope=('All125marker slots and70closed candidate species-tree views across526taxon entries retained; '
        '673complete-identical-protein predictor overlap cells freshly checked. Common original alignment '
        'coordinates and joint confidence/PAE-qualified observation mask, per-taxon min50/30%original-column '
        'eligibility, identical AA data/unknown masks/taxa/topology for both predictors. All8750comparison '
        'dispositions and4,523,750original internal-branch slots accounted, including no-overlap and merged '
        'paths. Exact input reuse preserves each marker/view mapping. DendroPy split construction independently '
        'checked against Bio.Phylo pruning of every raw candidate tree represented in ready cases; full serialized '
        'FASTA/config/positions/mapping readback. Future seven native roles per unique ready input are unlaunched: '
        'LG+F+G4 AA andAF+G4,AF+F+G4,LLM+G4 for each predictor on fixed unrooted topology. This selected predictor '
        'control is not coverage of the full fungal lineage, accepted species/gene-copy orthology/root/model, '
        'physical displacement, calibrated error/evolutionary effect, ancestral posterior or aim completion. '
        'No native inference/GPU/new costs or changes/restarts to original scientific jobs.')
    resources=dict(checked_utc=datetime.now(timezone.utc).isoformat(),cpus=2,memory_gib=8,swap_gib=0,blas_threads=1,
        planning_output_gib=8,minimum_free_disk_gib=16,maximum_comparison_cases=8750,maximum_unique_inputs=8750,
        future_native_roles_upper_count=61250,actual_future_native_roles_pending_complete_input_closure=True,
        active_hours_per_stage=[.01,4],runtime_uncalibrated=True,finish_eta=None,gpu=False,new_cost_usd=0,
        available_memory_gib=psutil.virtual_memory().available/2**30,available_disk_gib=shutil.disk_usage('.').free/2**30,
        scope='Input preparation/readback/closure only, bounded cgroup memory and no swap. Output/time planning is conservative and uncalibrated; no native fitting launch authority or ETA.')
    plan=dict(schema=SCHEMA,pins=bindings,output=output,completion=completion,launch_inventory=inventory,
              resources=resources,scope=scope)
    create(plan_path,plan)
    producer=launch('matched-predictor-branch-inputs',[sys.executable,'scripts/prepare_matched_predictor_branch_inputs.py',
        '--plan',plan_path],[],plan_path,cpus=2,memory=8)
    reader=launch('matched-predictor-branch-inputs-readback',[sys.executable,'scripts/readback_matched_predictor_branch_inputs.py',
        '--plan',plan_path],[producer],plan_path,cpus=2,memory=8)
    completion_plan='metadata/matched_predictor_branch_inputs_completion_plan_20261003_v1.json'
    spec=dict(source_plan=plan_path,producer_receipt=output+'/receipt.json',independent_readback=output+'/readback.json',
        producer_status=PRODUCER_STATUS,reader_status=READER_STATUS,completed_status='complete_verified_full_matched_predictor_branch_inputs',
        summary_fields=SUMMARY_FIELDS,launches=[producer,reader],pins={path:sha(path) for path in [plan_path,producer,reader,
            'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},output=completion,scope=scope)
    create(completion_plan,spec)
    closer=launch('matched-predictor-branch-inputs-closure',[sys.executable,'scripts/close_full_triad_sequence_stage.py',
        '--plan',completion_plan],[producer,reader],completion_plan,cpus=2,memory=8)
    create(inventory,dict(status='launched_full_matched_predictor_branch_input_workflow',source_plan=plan_path,
        source_plan_sha256=sha(plan_path),launches=[producer,reader,closer],native_fits_launched=0,gpu=False,new_cost_usd=0))


if __name__=='__main__':main()
