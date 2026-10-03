#!/usr/bin/env python3
"""Launch the complete native replay only after frozen contracts and closure."""
import argparse
import json
from pathlib import Path
import shutil
import sys

import psutil

from independent_native_clade_mapping import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import launch, create
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text())
    for path,d in plan['pins'].items():assert sha(path)==d,path
    for key,status in [('grid_validation','passed_full_independent_native_clade_mapping_contracts')]:
        record=json.loads(Path(plan[key]).read_text());assert record['status']==status
        for path,d in record['source_hashes'].items():assert sha(path)==d,path
    grid=json.loads(Path(plan['grid_validation']).read_text())
    assert grid['full_chains']==1620 and grid['checked_chains']==1618 and grid['failed_chains']==2
    assert grid['intact_chains_in_unresolved_groups']==6 and grid['candidate_frame_mappings']==653672
    assert len(grid['rehashed_false_exports_rejected'])==12
    assert len(grid['source_mapping_alterations_rejected'])==4
    assert grid['unchanged_interrupted_checkpoints_reused'] is True
    assert grid['completed_producer_and_alternate_reader_restart_refused'] is True
    assert psutil.virtual_memory().available>=32*2**30
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    assert plan['resources']['cpu_per_stage']==2 and plan['resources']['memory_gib']==32
    assert plan['resources']['swap_gib']==0 and plan['resources']['blas_threads']==1
    root=Path(plan['output']);assert not root.exists()
    dependencies=[plan['inventory_launch'],plan['diagnostic_closure_launch']]
    for path in dependencies:
        record=json.loads(Path(path).read_text());record['launch']=path
        assert fingerprint(record) is None;journal_terminal(record)
    producer=launch('independent-native-clade-mapping','metadata/independent_native_clade_mapping_wait_plan_20261002.json',
        [sys.executable,'scripts/prepare_independent_native_clade_mapping.py','--plan',str(a.plan)],dependencies,str(a.plan))
    reader=launch('independent-native-clade-mapping-readback','metadata/independent_native_clade_mapping_readback_wait_plan_20261002.json',
        [sys.executable,'scripts/readback_independent_native_clade_mapping.py','--plan',str(a.plan),'--output',str(root/'readback.json')],
        [producer],str(a.plan))
    completion=dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status='complete_full_independent_native_clade_mapping_pending_readback',
        reader_status='passed_full_independent_native_clade_mapping_serialized_readback',
        completed_status='complete_verified_full_independent_native_clade_mapping',summary_fields=SUMMARY_FIELDS,
        launches=[producer,reader],pins={p:sha(p) for p in [str(a.plan),producer,reader,'scripts/close_full_triad_sequence_stage.py',
        'scripts/record_completed_process_handoffs_v2.py']},output=plan['completion'],scope=plan['scope'])
    create(plan['completion_plan'],completion)
    closer=launch('independent-native-clade-mapping-closure','metadata/independent_native_clade_mapping_closure_wait_plan_20261002.json',
        [sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',plan['completion_plan']],[producer,reader],plan['completion_plan'])
    create(plan['launch_inventory'],dict(status='queued_full_independent_native_clade_mapping',source_plan=str(a.plan),source_plan_sha256=sha(a.plan),
        launches=[producer,reader,closer],native_sampling_restarted=False,gpu=False,new_cost_usd=0,scope=plan['scope']))


if __name__=='__main__':main()
