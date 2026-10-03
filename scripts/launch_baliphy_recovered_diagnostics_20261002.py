#!/usr/bin/env python3
"""Launch the complete closed recovered-horizon diagnostic publication."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from launch_full_triad_sequence_geometry_followup import launch
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan',type=Path,required=True); a=p.parse_args()
    plan=json.loads(a.plan.read_text())
    for path,digest in plan['pins'].items(): assert sha(path)==digest,path
    assert psutil.virtual_memory().available >= 32*2**30
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib']*2**30
    from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
    dependency=json.loads(Path(plan['diagnostic_closure_launch']).read_text()); dependency['launch']=plan['diagnostic_closure_launch']
    assert fingerprint(dependency) is None
    journal_terminal(dependency)
    assert not Path(plan['table']).exists() and not Path(plan['completion']).exists()
    launch('baliphy-recovered-diagnostics-publication', 'metadata/baliphy_recovered_diagnostics_publication_wait_plan_20261002.json',
        [sys.executable,'scripts/publish_baliphy_recovered_diagnostics_20261002.py','--plan',str(a.plan)],
        [plan['diagnostic_closure_launch']], str(a.plan))


if __name__ == '__main__': main()
