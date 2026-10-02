#!/usr/bin/env python3
"""Queue original-context directions behind the original joint-direction closure."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from full_triad_context_joint_direction_sources import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import launch,create
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text())
    for path,digest in plan['pins'].items():assert sha(path)==digest,path
    dependency=plan['joint_closure_launch'];r=json.loads(Path(dependency).read_text());proc=psutil.Process(r['pid'])
    assert proc.create_time()==r['created'] and proc.cmdline()==r['cmdline']
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30 and psutil.virtual_memory().available>=32*2**30
    f=json.loads(Path(plan['fixture_validation']).read_text());assert f['status']=='passed_full_triad_context_joint_direction_software_contracts' and f['rejected_false_exports']==23 and sha(f['receipt'])==f['receipt_sha256']
    name='full-triad-context-joint-directions';prefix='metadata/full_triad_context_joint_directions';root=Path(plan['output'])
    producer=launch(name,prefix+'_wait_plan_20261002.json',[sys.executable,'scripts/project_full_triad_context_joint_directions.py','--plan',str(a.plan)],[dependency],str(a.plan))
    reader=launch(name+'-readback',prefix+'_readback_wait_plan_20261002.json',[sys.executable,'scripts/readback_full_triad_context_joint_directions.py','--plan',str(a.plan),'--output',str(root/'readback.json')],[producer],str(a.plan))
    cp=prefix+'_completion_plan_20261002.json'
    create(cp,dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status='complete_full_triad_context_joint_directions_pending_independent_readback',reader_status='passed_full_triad_context_joint_directions_original_gate_sql_readback',
        completed_status='complete_verified_full_triad_context_joint_directions',summary_fields=SUMMARY_FIELDS,launches=[producer,reader],
        pins={str(p):sha(p) for p in [a.plan,producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=prefix+'_completed_20261002.json',scope=plan['scope']))
    closer=launch(name+'-closure',prefix+'_closure_wait_plan_20261002.json',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],[producer,reader],cp)
    create(prefix+'_launches_20261002.json',dict(status='queued_full_triad_context_joint_directions',source_plan=str(a.plan),source_plan_sha256=sha(a.plan),launches=[producer,reader,closer],scope=plan['scope']))


if __name__=='__main__':main()
