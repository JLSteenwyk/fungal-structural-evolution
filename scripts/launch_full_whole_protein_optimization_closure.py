#!/usr/bin/env python3
"""Queue exhaustive optimization closure behind the exact original full reader."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from launch_full_triad_sequence_geometry_followup import launch
from screen_duplication_alignment_reuse import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text())
    for p,h in plan['pins'].items():assert sha(p)==h,p
    f=json.loads(Path(plan['fixture_validation']).read_text());assert f['status']=='passed_full_whole_protein_optimization_closure_contracts' and f['rejected_false_evidence_cases']==11
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30 and psutil.virtual_memory().available>=32*2**30
    dependency=plan['launches'][-1];r=json.loads(Path(dependency).read_text());p=psutil.Process(r['pid'])
    assert p.create_time()==r['created'] and p.cmdline()==r['cmdline'] and p.status()!=psutil.STATUS_ZOMBIE
    launch('full-whole-protein-optimization-closure','metadata/full_whole_protein_optimization_closure_wait_plan_20261002.json',
        [sys.executable,'scripts/close_full_whole_protein_optimization.py','--plan',str(args.plan)],[dependency],str(args.plan))


if __name__=='__main__':main()
