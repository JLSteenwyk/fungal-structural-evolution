#!/usr/bin/env python3
"""Launch full logical-case indexing and independent readback of closed matching."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from full_matching_case_sources import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import launch,create
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text())
    for path,digest in plan['pins'].items():assert sha(path)==digest,path
    source=json.loads(Path(plan['matching_completion']).read_text())
    assert source['status']=='complete_verified_full_fixed_matched_coverage_attrition' and source['exact_process_journals_checked']==2
    assert sha(source['full_hash_archive'])==source['full_hash_archive_sha256']
    f=json.loads(Path(plan['fixture_validation']).read_text())
    assert f['status']=='passed_full_matching_case_index_software_contracts' and len(f['rejected_rehashed_exports'])==23
    assert f['distinct_genes_sharing_models_preserved'] and f['completed_restart_refused'] and f['interrupted_full_replay_passed']
    for path,digest in f['source_hashes'].items():assert sha(path)==digest,path
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30 and psutil.virtual_memory().available>=32*2**30
    root=Path(plan['output']);prefix='metadata/full_matching_case_index';name='full-matching-case-index'
    producer=launch(name,prefix+'_wait_plan_20261002.json',[sys.executable,'scripts/index_full_matching_cases.py','--plan',str(a.plan)],[],str(a.plan))
    reader=launch(name+'-readback',prefix+'_readback_wait_plan_20261002.json',[sys.executable,'scripts/readback_full_matching_cases.py','--plan',str(a.plan),'--output',str(root/'readback.json')],[producer],str(a.plan))
    cp=prefix+'_completion_plan_20261002.json'
    create(cp,dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),producer_status='complete_full_matching_case_index_pending_independent_readback',reader_status='passed_full_matching_case_index_sql_readback',completed_status='complete_verified_full_matching_logical_case_index',summary_fields=SUMMARY_FIELDS,launches=[producer,reader],pins={str(q):sha(q) for q in [a.plan,producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},output=prefix+'_completed_20261002.json',scope=plan['scope']))
    closer=launch(name+'-closure',prefix+'_closure_wait_plan_20261002.json',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],[producer,reader],cp)
    create(prefix+'_launches_20261002.json',dict(status='launched_full_matching_logical_case_index',source_plan=str(a.plan),source_plan_sha256=sha(a.plan),launches=[producer,reader,closer],scope=plan['scope']))


if __name__=='__main__':main()
