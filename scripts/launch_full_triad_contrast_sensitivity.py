#!/usr/bin/env python3
"""Launch frozen full-grid contrast summaries, raw-fit reader and journal closure."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from full_triad_contrast_sensitivity_sources import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import launch, create
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text())
    for path,digest in plan['pins'].items():assert sha(path)==digest,path
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    assert psutil.virtual_memory().available>=plan['resources']['memory_gib']*2**30
    fixture=json.loads(Path(plan['fixture_validation']).read_text());assert fixture['status']=='passed_full_triad_contrast_sensitivity_software_contracts'
    assert sha(fixture['receipt'])==fixture['receipt_sha256'] and fixture['rejected_false_exports']==16
    root=Path(plan['output']);name='full-triad-contrast-sensitivity';prefix='metadata/full_triad_contrast_sensitivity'
    producer=launch(name,prefix+'_wait_plan_20261002.json',
        [sys.executable,'scripts/summarize_full_triad_contrast_sensitivity.py','--plan',str(a.plan)],[],str(a.plan))
    reader=launch(name+'-readback',prefix+'_readback_wait_plan_20261002.json',
        [sys.executable,'scripts/readback_full_triad_contrast_sensitivity.py','--plan',str(a.plan),'--output',str(root/'readback.json')],[producer],str(a.plan))
    cp=prefix+'_completion_plan_20261002.json'
    create(cp,dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status='complete_full_triad_contrast_sensitivity_pending_independent_readback',reader_status='passed_full_triad_contrast_sensitivity_raw_fit_sql_readback',
        completed_status='complete_verified_full_triad_contrast_sensitivity',summary_fields=SUMMARY_FIELDS,launches=[producer,reader],
        pins={str(p):sha(p) for p in [a.plan,producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=prefix+'_completed_20261002.json',scope=plan['scope']))
    closer=launch(name+'-closure',prefix+'_closure_wait_plan_20261002.json',
        [sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],[producer,reader],cp)
    create(prefix+'_launches_20261002.json',dict(status='launched_full_triad_contrast_sensitivity',source_plan=str(a.plan),source_plan_sha256=sha(a.plan),
        launches=[producer,reader,closer],scope=plan['scope']))


if __name__=='__main__':main()
