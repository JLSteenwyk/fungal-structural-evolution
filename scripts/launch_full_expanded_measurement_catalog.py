#!/usr/bin/env python3
"""Launch the entire closed target/background catalog and independent reader."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from full_expanded_measurement_catalog_sources import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import launch,create
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text())
    for path,d in plan['pins'].items():assert sha(path)==d,path
    c=json.loads(Path(plan['matching_completion']).read_text());assert c['status']=='complete_verified_full_fixed_matched_coverage_attrition' and sha(c['full_hash_archive'])==c['full_hash_archive_sha256']
    f=json.loads(Path(plan['fixture_validation']).read_text());assert f['status']=='passed_full_expanded_measurement_catalog_software_contracts' and len(f['rejected_rehashed_exports'])==20 and f['completed_restart_refused'] and f['interrupted_full_replay_passed']
    for path,d in f['source_hashes'].items():assert sha(path)==d,path
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30 and psutil.virtual_memory().available>=32*2**30
    root=Path(plan['output']);prefix='metadata/full_expanded_measurement_catalog';name='full-expanded-measurement-catalog'
    producer=launch(name,prefix+'_wait_plan_20261002.json',[sys.executable,'scripts/export_full_expanded_measurement_catalog.py','--plan',str(a.plan)],[],str(a.plan))
    reader=launch(name+'-readback',prefix+'_readback_wait_plan_20261002.json',[sys.executable,'scripts/readback_full_expanded_measurement_catalog.py','--plan',str(a.plan),'--output',str(root/'readback.json')],[producer],str(a.plan))
    cp=prefix+'_completion_plan_20261002.json'
    create(cp,dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),producer_status='complete_full_expanded_measurement_catalog_pending_independent_readback',reader_status='passed_full_expanded_measurement_catalog_source_and_numeric_readback',completed_status='complete_verified_full_expanded_directed_measurement_catalog',summary_fields=SUMMARY_FIELDS,launches=[producer,reader],pins={str(q):sha(q) for q in [a.plan,producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},output=prefix+'_completed_20261002.json',scope=plan['scope']))
    closer=launch(name+'-closure',prefix+'_closure_wait_plan_20261002.json',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],[producer,reader],cp)
    create(prefix+'_launches_20261002.json',dict(status='launched_full_expanded_directed_measurement_catalog',source_plan=str(a.plan),source_plan_sha256=sha(a.plan),launches=[producer,reader,closer],scope=plan['scope']))


if __name__=='__main__':main()
