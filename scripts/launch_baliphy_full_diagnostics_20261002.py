#!/usr/bin/env python3
"""Launch the complete closed first-horizon diagnostic publication."""
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
    assert not Path(plan['table']).exists() and not Path(plan['completion']).exists()
    launch('baliphy-full-diagnostics-publication', 'metadata/baliphy_full_diagnostics_publication_wait_plan_20261002.json',
        [sys.executable,'scripts/publish_baliphy_full_diagnostics_20261002.py','--plan',str(a.plan)],
        [plan['postprocessing_accounting_launch']], str(a.plan))


if __name__ == '__main__': main()
