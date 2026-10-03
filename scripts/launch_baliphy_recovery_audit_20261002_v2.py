#!/usr/bin/env python3
"""Launch read-only recovery accounting; never restart a native chain."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from launch_full_triad_sequence_geometry_followup import create, launch
from record_project_runtime_checkpoint_v4 import journal_terminal
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text())
    for path,digest in plan['pins'].items():assert sha(path)==digest,path
    check=json.loads(Path(plan['validation']).read_text())
    assert check['status']=='passed_full_baliphy_recovery_overlay_contracts'
    assert len(check['malformed_overlays_rejected'])==18
    for path,digest in check['source_hashes'].items():assert sha(path)==digest,path
    assert not Path(plan['output']).exists() and not Path(plan['completion']).exists()
    assert not Path(plan['closure_output']).exists() and not Path(plan['closure_plan']).exists()
    assert psutil.virtual_memory().available>=32*2**30
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    original=json.loads(Path(plan['recovery_launch']).read_text());original['launch']=plan['recovery_launch']
    journal_terminal(original)
    reader=launch('baliphy-memory-recovery-readback-v2',
        'metadata/baliphy_memory_recovery_readback_wait_plan_20261002_v2.json',
        [sys.executable,'scripts/audit_baliphy_memory_recovery_20261002_v2.py','--plan',str(args.plan)],
        [plan['recovery_launch']],str(args.plan))
    pins={str(p):sha(p) for p in [args.plan,reader,plan['recovery_launch'],
        'scripts/close_baliphy_memory_recovery_20261002.py','scripts/record_completed_process_handoffs_v2.py',
        'scripts/reference_measurement_union_sources.py','scripts/run_ortholog_pair_guide_comparison.py']}
    create(plan['closure_plan'],dict(audit_plan=str(args.plan),launches=[plan['recovery_launch'],reader],
        pins=pins,output=plan['closure_output'],completion=plan['completion'],scope=plan['scope']))
    closer=launch('baliphy-memory-recovery-accounting-v2',
        'metadata/baliphy_memory_recovery_accounting_wait_plan_20261002_v2.json',
        [sys.executable,'scripts/close_baliphy_memory_recovery_20261002.py','--plan',plan['closure_plan']],
        [reader],plan['closure_plan'])
    print(json.dumps(dict(status='queued_full_original_grid_recovery_readback_and_accounting',
        reader_launch=reader,closure_launch=closer,native_restarted=False)),flush=True)


if __name__=='__main__':main()
