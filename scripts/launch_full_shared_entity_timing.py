#!/usr/bin/env python3
"""Queue whole-grid technical timing after the exact original qualification."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from launch_full_triad_sequence_geometry_followup import create,launch
from record_project_runtime_checkpoint_v4 import fingerprint,journal_terminal
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text())
    for path,digest in plan['pins'].items():assert sha(path)==digest,path
    checked=json.loads(Path(plan['validation']).read_text())
    assert checked['status']=='passed_full_shared_entity_timing_grid_checkpoint_and_numerical_probe_contracts'
    assert checked['synthetic_cases']==7200 and checked['synthetic_cohorts']==6
    assert len(checked['rehashed_census_alterations_rejected'])==len(checked['malformed_numeric_probe_exports_rejected'])==6
    for path,digest in checked['source_hashes'].items():assert sha(path)==digest,path
    assert sha(plan['fit_plan'])==plan['fit_plan_sha256']
    fit=json.loads(Path(plan['fit_plan']).read_text());assert fit['launch_state']=='not_launched_or_queued'
    assert not Path(fit['output']).exists() and not Path(plan['output']).exists() and not Path(plan['completion']).exists()
    assert not Path(plan['closure_plan']).exists()
    assert psutil.virtual_memory().available>=32*2**30
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    original=json.loads(Path(plan['qualification_closure_launch']).read_text());original['launch']=plan['qualification_closure_launch']
    assert sha(original['plan'])==original['plan_sha256']
    if fingerprint(original) is None:journal_terminal(original)
    producer=launch('full-shared-entity-timing','metadata/full_shared_entity_timing_wait_plan_20261002.json',
        [sys.executable,'scripts/prepare_full_shared_entity_timing.py','--plan',str(a.plan)],
        [plan['qualification_closure_launch']],str(a.plan))
    root=Path(plan['output'])
    reader=launch('full-shared-entity-timing-readback','metadata/full_shared_entity_timing_readback_wait_plan_20261002.json',
        [sys.executable,'scripts/readback_full_shared_entity_timing.py','--plan',str(a.plan),'--output',str(root/'readback.json')],
        [producer],str(a.plan))
    spec=dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status='complete_full_shared_entity_timing_pending_accounting_readback',
        reader_status='passed_full_shared_entity_timing_census_and_planning_accounting_readback',
        completed_status='complete_verified_full_shared_entity_timing_accounting',summary_fields=plan['summary_fields'],
        launches=[producer,reader],pins={str(p):sha(p) for p in [a.plan,producer,reader,
            'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'],scope=plan['scope'])
    create(plan['closure_plan'],spec)
    closer=launch('full-shared-entity-timing-closure','metadata/full_shared_entity_timing_closure_wait_plan_20261002.json',
        [sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',plan['closure_plan']],
        [producer,reader],plan['closure_plan'])
    create(plan['launch_inventory'],dict(status='queued_full_shared_entity_qualified_input_timing',
        source_plan=str(a.plan),source_plan_sha256=sha(a.plan),launches=[producer,reader,closer],
        production_fitting_launched=False,gpu=False,new_cost_usd=0,scope=plan['scope']))


if __name__=='__main__':main()
