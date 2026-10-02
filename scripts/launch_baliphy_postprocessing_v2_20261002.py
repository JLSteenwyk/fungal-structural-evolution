#!/usr/bin/env python3
"""Launch the corrected full accounting after saved-schema regression checks."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from launch_full_triad_sequence_geometry_followup import launch
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args(); plan = json.loads(a.plan.read_text())
    for path, digest in plan['pins'].items(): assert sha(path) == digest, path
    v = json.loads(Path(plan['coordinate_validation']).read_text())
    assert v['status'] == 'passed_full_saved_quartet_coordinate_schema_contract'
    assert v['complete_quartets'] == 402 and v['chain_coordinate_metadata_checked'] == 1608
    assert v['cutoff_summaries_checked'] == v['original_bad_count_reproduced'] == 804
    assert len(v['malformed_metadata_rejected']) == 9
    for path, digest in v['script_pins'].items(): assert sha(path) == digest, path
    assert not Path(plan['output']).exists() and not Path(plan['completion']).exists()
    assert psutil.virtual_memory().available >= 32 * 2**30
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    launch('baliphy-full-postprocessing-accounting-v2',
        'metadata/baliphy_full_postprocessing_accounting_wait_plan_20261002_v2.json',
        [sys.executable, 'scripts/close_baliphy_postprocessing_20261002_v2.py', '--plan', str(a.plan)], [], str(a.plan))


if __name__ == '__main__': main()
