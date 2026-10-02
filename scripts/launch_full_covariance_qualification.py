#!/usr/bin/env python3
"""Queue complete uniform covariance audits behind the original design closure."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from full_covariance_qualification_sources import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import launch,create
from record_project_runtime_checkpoint_v4 import fingerprint,journal_terminal
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text())
    for path,d in plan['pins'].items():assert sha(path)==d,path
    fixture=json.loads(Path(plan['fixture_validation']).read_text());primitive=json.loads(Path(plan['primitive_validation']).read_text())
    assert fixture['status']=='passed_full_uniform_covariance_qualification_contracts' and len(fixture['rejected_rehashed_exports'])==17
    assert primitive['status']=='passed_reusable_covariance_kernel_context_contracts' and primitive['independent_dense_and_latent_cases']==48
    for report in [fixture,primitive]:
        for path,d in report['source_hashes'].items():assert sha(path)==d,path
    for key in ['full_interrupt_replay_passed','completed_restart_refused','source_and_journal_fixtures_synthetic',
                'all_loading_modes_trees_empty_and_review_cases_retained','independent_seven_basis_qualified_path_passed']:assert fixture[key]
    for key,status in [('operator_completion','complete_verified_full_entity_operator_bank'),('inputs_completion','complete_verified_full_expanded_model_inputs')]:
        r=json.loads(Path(plan[key]).read_text());assert r['status']==status and r['exact_process_journals_checked']==2
        assert sha(r['full_hash_archive'])==r['full_hash_archive_sha256']
    for path in plan['dependencies']:
        r=json.loads(Path(path).read_text());r['launch']=path
        assert sha(r['plan'])==r['plan_sha256']
        if fingerprint(r) is None:journal_terminal(r)
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    assert psutil.virtual_memory().available>=32*2**30
    root=Path(plan['output']);assert root.parent.is_dir() and not root.exists()
    prefix='metadata/full_uniform_covariance_qualification';name='full-uniform-covariance-qualification'
    producer=launch(name,prefix+'_wait_plan_20261002.json',[sys.executable,'scripts/prepare_full_covariance_qualification.py','--plan',str(a.plan)],plan['dependencies'],str(a.plan))
    reader=launch(name+'-readback',prefix+'_readback_wait_plan_20261002.json',[sys.executable,'scripts/readback_full_covariance_qualification.py','--plan',str(a.plan),
        '--output',str(root/'readback.json')],[producer],str(a.plan))
    cp=prefix+'_completion_plan_20261002.json'
    create(cp,dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status='complete_full_uniform_covariance_qualification_pending_independent_readback',reader_status='passed_full_uniform_covariance_latent_and_sql_readback',
        completed_status='complete_verified_full_uniform_covariance_qualification',summary_fields=SUMMARY_FIELDS,launches=[producer,reader],
        pins={str(x):sha(x) for x in [a.plan,producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=prefix+'_completed_20261002.json',scope=plan['scope']))
    closure=launch(name+'-closure',prefix+'_closure_wait_plan_20261002.json',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],[producer,reader],cp)
    create(prefix+'_launches_20261002.json',dict(status='queued_full_uniform_covariance_qualification',source_plan=str(a.plan),source_plan_sha256=sha(a.plan),
        launches=[producer,reader,closure],scope=plan['scope']))


if __name__=='__main__':main()
