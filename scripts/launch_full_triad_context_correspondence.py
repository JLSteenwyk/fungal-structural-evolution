#!/usr/bin/env python3
"""Queue the complete original-context correspondence join after original comparison closure."""
import argparse
import json
from pathlib import Path
import sys
import psutil
from launch_full_triad_sequence_geometry_followup import launch,create
from full_triad_context_correspondence_sources import SUMMARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workflow',type=Path,required=True);args=p.parse_args()
    workflow=json.loads(args.workflow.read_text())
    for path,digest in workflow['pins'].items():assert sha(path)==digest
    dependency=workflow['comparison_closure_launch'];original=json.loads(Path(dependency).read_text())
    proc=psutil.Process(original['pid']);assert proc.create_time()==original['created'] and proc.cmdline()==original['cmdline']
    source=workflow['source_plan'];config=json.loads(Path(source).read_text());root=Path(config['output']);prefix='metadata/full_triad_context_correspondence';name='full-triad-context-correspondence'
    producer=launch(name,prefix+'_wait_plan_20261002.json',
        [sys.executable,'scripts/project_full_triad_context_correspondence.py','--plan',source],[dependency],source)
    reader=launch(name+'-readback',prefix+'_readback_wait_plan_20261002.json',
        [sys.executable,'scripts/readback_full_triad_context_correspondence.py','--plan',source,'--output',str(root/'readback.json')],[producer],source)
    close_plan=prefix+'_completion_plan_20261002.json'
    create(close_plan,dict(source_plan=source,producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status='complete_full_triad_context_correspondence_pending_independent_readback',reader_status='passed_full_triad_context_correspondence_sql_readback',
        completed_status='complete_verified_full_triad_context_correspondence',summary_fields=SUMMARY_FIELDS,launches=[producer,reader],
        pins={str(p):sha(p) for p in [source,producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output='metadata/full_triad_context_correspondence_completed_20261002.json',scope=config['scope']))
    closer=launch(name+'-closure',prefix+'_closure_wait_plan_20261002.json',
        [sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',close_plan],[producer,reader],close_plan)
    create(workflow['launch_inventory'],dict(status='queued_full_original_context_correspondence_projection',workflow_sha256=sha(args.workflow),launches=[producer,reader,closer],scope=config['scope']))


if __name__=='__main__':main()
