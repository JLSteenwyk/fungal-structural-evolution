#!/usr/bin/env python3
"""Queue full short joint output qualification behind four original closures."""
import argparse
import json
from pathlib import Path
import shutil
import sys

import psutil

from ancestral_chain_attempt import sha
from baliphy_joint_sampler_qualification_v3 import SUMMARY_FIELDS
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from record_project_runtime_checkpoint_v4 import fingerprint,journal_terminal
from run_baliphy_reference_preflight import verify
from run_baliphy_joint_sampler_qualification_v3 import PRODUCER_STATUS,READER_STATUS,COMPLETED_STATUS


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());verify(plan)
    gate=json.loads(Path(plan['software_validation']).read_text())
    assert gate['status']=='passed_full_joint_sampler_qualification_v3_software_contracts'
    verify({'pins':gate['source_hashes']})
    assert gate['full_roles']==1620 and gate['full_quartets']==405
    assert {'missing_frame_export','missing_export_directory'} <= set(gate['serialization_cases_rejected'])
    assert plan['resources']['iterations']==20 and plan['resources']['posterior_qualified'] is False
    assert len(plan['dependencies'])==4
    assert psutil.virtual_memory().available>=200*2**30 and shutil.disk_usage('.').free>=256*2**30
    assert not Path(plan['output']).exists()
    for path in plan['dependencies']:
        record=json.loads(Path(path).read_text());record['launch']=path
        assert sha(record['plan'])==record['plan_sha256']
        if fingerprint(record) is None:journal_terminal(record)
    producer=launch('baliphy-joint-sampler-qualification-v3',
        [sys.executable,'scripts/run_baliphy_joint_sampler_qualification_v3.py','--plan',str(a.plan)],
        plan['dependencies'],str(a.plan),cpus=16,memory=200)
    reader=launch('baliphy-joint-sampler-qualification-v3-readback',
        [sys.executable,'scripts/run_baliphy_joint_sampler_qualification_v3.py','--plan',str(a.plan),'--reader'],
        [producer],str(a.plan),cpus=2,memory=32)
    root=Path(plan['output'])
    completion=dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),
        independent_readback=str(root/'readback.json'),producer_status=PRODUCER_STATUS,reader_status=READER_STATUS,
        completed_status=COMPLETED_STATUS,summary_fields=SUMMARY_FIELDS+['reservation_audit'],launches=[producer,reader],
        pins={p:sha(p) for p in [str(a.plan),producer,reader,'scripts/close_full_triad_sequence_stage.py',
            'scripts/record_completed_process_handoffs_v2.py']},output=plan['completion'],scope=plan['scope'])
    create(plan['completion_plan'],completion)
    closer=launch('baliphy-joint-sampler-qualification-v3-closure',
        [sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',plan['completion_plan']],
        [producer,reader],plan['completion_plan'],cpus=2,memory=32)
    create(plan['launch_inventory'],dict(status='queued_full_joint_short_sampler_qualification_v3',
        source_plan=str(a.plan),source_plan_sha256=sha(a.plan),launches=[producer,reader,closer],
        original_prerequisite_closures=plan['dependencies'],expected_native_roles=1620,iterations_per_role=20,
        longer_posterior_horizon_launched=False,posterior_qualified=False,gpu=False,new_cost_usd=0))


if __name__=='__main__':main()
