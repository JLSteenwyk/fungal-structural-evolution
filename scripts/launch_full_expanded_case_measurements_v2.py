#!/usr/bin/env python3
"""Queue the entire matched measurement join after both original source closures."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from full_expanded_case_measurement_sources import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import launch,create
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text())
    for path,d in plan['pins'].items():assert sha(path)==d,path
    dependencies=[plan[k+'_closure_launch'] for k in ['case_index','catalog']]
    for path in dependencies:
        r=json.loads(Path(path).read_text());proc=psutil.Process(r['pid']);assert proc.create_time()==r['created'] and proc.cmdline()==r['cmdline'] and proc.status()!=psutil.STATUS_ZOMBIE
    f=json.loads(Path(plan['fixture_validation']).read_text());assert f['status']=='passed_full_expanded_case_measurement_join_software_contracts' and len(f['rejected_rehashed_exports'])==19 and f['independent_decimal_arithmetic_passed'] and f['interrupted_full_replay_passed'] and f['completed_restart_refused']
    for path,d in f['source_hashes'].items():assert sha(path)==d,path
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30 and psutil.virtual_memory().available>=32*2**30
    root=Path(plan['output']);prefix='metadata/full_expanded_case_measurements_v2';name='full-expanded-case-measurements-v2'
    producer=launch(name,prefix+'_wait_plan_20261002.json',[sys.executable,'scripts/join_full_expanded_case_measurements.py','--plan',str(a.plan)],dependencies,str(a.plan))
    reader=launch(name+'-readback',prefix+'_readback_wait_plan_20261002.json',[sys.executable,'scripts/readback_full_expanded_case_measurements.py','--plan',str(a.plan),'--output',str(root/'readback.json')],[producer],str(a.plan))
    cp=prefix+'_completion_plan_20261002.json'
    create(cp,dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),producer_status='complete_full_expanded_case_measurement_join_pending_independent_readback',reader_status='passed_full_expanded_case_measurement_join_decimal_readback',completed_status='complete_verified_full_expanded_matched_case_measurements',summary_fields=SUMMARY_FIELDS,launches=[producer,reader],pins={str(q):sha(q) for q in [a.plan,producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},output=prefix+'_completed_20261002.json',scope=plan['scope']))
    closer=launch(name+'-closure',prefix+'_closure_wait_plan_20261002.json',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],[producer,reader],cp)
    create(prefix+'_launches_20261002.json',dict(status='queued_full_expanded_matched_case_measurements',source_plan=str(a.plan),source_plan_sha256=sha(a.plan),launches=[producer,reader,closer],scope=plan['scope']))


if __name__=='__main__':main()
