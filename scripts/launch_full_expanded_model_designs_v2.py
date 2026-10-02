#!/usr/bin/env python3
"""Queue full rank/recipe inventory behind the original complete input closure."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from full_expanded_model_design_sources import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import launch, create
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text())
    for path,d in plan['pins'].items():assert sha(path)==d,path
    fixture=json.loads(Path(plan['fixture_validation']).read_text())
    assert fixture['status']=='passed_full_expanded_model_design_v2_software_contracts'
    assert len(fixture['rejected_rehashed_exports'])==23
    for k in ['full_interrupt_replay_passed','completed_restart_refused','extreme_scale_qr_svd_check_passed','full_grid_retained_and_exact_cohorts_shared','source_and_journal_fixtures_synthetic','rank_boundary_condition_readback_passed']:assert fixture[k]
    for path,d in fixture['source_hashes'].items():assert sha(path)==d,path
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    assert psutil.virtual_memory().available>=32*2**30
    dependency=plan['inputs_closure_launch'];record=json.loads(Path(dependency).read_text());record['launch']=dependency
    if fingerprint(record) is None:
        journal_terminal(record)
        assert json.loads(Path(plan['inputs_completion']).read_text())['status']=='complete_verified_full_expanded_model_inputs'
    root=Path(plan['output']);assert root.parent.is_dir() and not root.exists()
    prefix='metadata/full_expanded_model_designs_v2';name='full-expanded-model-designs-v2'
    producer=launch(name,prefix+'_wait_plan_20261002.json',
        [sys.executable,'scripts/prepare_full_expanded_model_designs.py','--plan',str(a.plan)],[dependency],str(a.plan))
    reader=launch(name+'-readback',prefix+'_readback_wait_plan_20261002.json',
        [sys.executable,'scripts/readback_full_expanded_model_designs_v2.py','--plan',str(a.plan),'--output',str(root/'readback.json')],[producer],str(a.plan))
    cp=prefix+'_completion_plan_20261002.json'
    create(cp,dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status='complete_full_expanded_model_designs_pending_independent_readback',reader_status='passed_full_expanded_model_design_sql_qr_readback',
        completed_status='complete_verified_full_expanded_model_designs',summary_fields=SUMMARY_FIELDS,launches=[producer,reader],
        pins={str(x):sha(x) for x in [a.plan,producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=prefix+'_completed_20261002.json',scope=plan['scope']))
    closure=launch(name+'-closure',prefix+'_closure_wait_plan_20261002.json',
        [sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],[producer,reader],cp)
    create(prefix+'_launches_20261002.json',dict(status='queued_full_expanded_model_designs',source_plan=str(a.plan),source_plan_sha256=sha(a.plan),launches=[producer,reader,closure],scope=plan['scope']))


if __name__=='__main__':main()
