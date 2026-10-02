#!/usr/bin/env python3
"""Launch the complete closed-source entity bank, independent readback and closure."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from full_entity_operator_sources import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import launch,create
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text())
    for path,d in plan['pins'].items():assert sha(path)==d,path
    c=json.loads(Path(plan['covariance_completion']).read_text());assert c['status']=='complete_verified_full_expanded_covariance' and c['exact_process_journals_checked']==2
    assert sha(c['full_hash_archive'])==c['full_hash_archive_sha256']
    fixture=json.loads(Path(plan['fixture_validation']).read_text());primitive=json.loads(Path(plan['primitive_validation']).read_text())
    assert fixture['status']=='passed_full_entity_operator_bank_software_contracts' and len(fixture['rejected_rehashed_exports'])==15
    assert primitive['status']=='passed_shared_entity_block_low_rank_covariance_contracts' and primitive['dense_cases']==36 and primitive['high_precision_signed_shared_case_passed']
    for item in [fixture,primitive]:
        for path,d in item['source_hashes'].items():assert sha(path)==d,path
    for k in ['full_interrupt_replay_passed','completed_restart_refused','source_and_journal_fixtures_synthetic','all_five_trees_both_loadings_preserved','distinct_genes_shared_models_and_cancelling_family_preserved']:assert fixture[k]
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30 and psutil.virtual_memory().available>=32*2**30
    root=Path(plan['output']);assert root.parent.is_dir() and not root.exists()
    prefix='metadata/full_entity_operator_bank';name='full-entity-operator-bank'
    producer=launch(name,prefix+'_wait_plan_20261002.json',[sys.executable,'scripts/prepare_full_entity_operators.py','--plan',str(a.plan)],[],str(a.plan))
    reader=launch(name+'-readback',prefix+'_readback_wait_plan_20261002.json',[sys.executable,'scripts/readback_full_entity_operators.py','--plan',str(a.plan),'--output',str(root/'readback.json')],[producer],str(a.plan))
    cp=prefix+'_completion_plan_20261002.json'
    create(cp,dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status='complete_full_entity_operator_bank_pending_independent_readback',reader_status='passed_full_entity_operator_raw_loading_and_gram_readback',completed_status='complete_verified_full_entity_operator_bank',
        summary_fields=SUMMARY_FIELDS,launches=[producer,reader],pins={str(x):sha(x) for x in [a.plan,producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=prefix+'_completed_20261002.json',scope=plan['scope']))
    closure=launch(name+'-closure',prefix+'_closure_wait_plan_20261002.json',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],[producer,reader],cp)
    create(prefix+'_launches_20261002.json',dict(status='launched_full_entity_operator_bank',source_plan=str(a.plan),source_plan_sha256=sha(a.plan),launches=[producer,reader,closure],scope=plan['scope']))


if __name__=='__main__':main()
