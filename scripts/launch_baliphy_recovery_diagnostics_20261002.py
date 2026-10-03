#!/usr/bin/env python3
"""Launch full-grid recovery postprocessing, readback and provenance closure."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import psutil
from launch_full_triad_sequence_geometry_followup import create,launch
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text())
    for path,digest in plan['pins'].items():assert sha(path)==digest,path
    check=json.loads(Path(plan['validation']).read_text())
    assert check['status']=='passed_full_baliphy_recovery_diagnostic_accounting_and_input_contracts'
    assert check['quartet_routes']==dict(reuse_closed_original_unchanged_quartet=402,compute_complete_recovery_quartet=1,unresolved_failed_native_chain_retained=2)
    assert len(check['malformed_new_inputs_rejected'])==6 and check['foreign_unchanged_quartet_report_rejected'] is True
    for path,digest in check['source_hashes'].items():assert sha(path)==digest,path
    completed=json.loads(Path(plan['recovery_completion']).read_text())
    assert completed['status']=='complete_verified_full_baliphy_memory_recovery_accounting' and completed['complete_quartets']==403
    assert psutil.virtual_memory().available>=32*2**30
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    assert not Path(plan['output']).exists() and not Path(plan['completion']).exists()
    assert not Path(plan['closure_plan']).exists()
    producer=launch('baliphy-recovery-full-diagnostics','metadata/baliphy_recovery_full_diagnostics_wait_plan_20261002.json',
        [sys.executable,'scripts/prepare_baliphy_recovery_diagnostics_20261002.py','--plan',str(a.plan)],
        [plan['recovery_closure_launch']],str(a.plan))
    root=Path(plan['output'])
    reader=launch('baliphy-recovery-full-diagnostics-readback','metadata/baliphy_recovery_full_diagnostics_readback_wait_plan_20261002.json',
        [sys.executable,'scripts/readback_baliphy_recovery_diagnostics_20261002.py','--plan',str(a.plan),
         '--output',str(root/'readback.json')],[producer],str(a.plan))
    spec=dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status='complete_full_baliphy_recovery_diagnostics_pending_accounting_readback',
        reader_status='passed_full_baliphy_recovery_diagnostic_accounting_and_new_quartet_inputs',
        completed_status='complete_verified_full_baliphy_recovery_diagnostics',summary_fields=plan['summary_fields'],
        launches=[producer,reader],pins={str(p):sha(p) for p in [a.plan,producer,reader,
            'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'],scope=plan['scope'])
    create(plan['closure_plan'],spec)
    closer=launch('baliphy-recovery-full-diagnostics-closure','metadata/baliphy_recovery_full_diagnostics_closure_wait_plan_20261002.json',
        [sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',plan['closure_plan']],
        [producer,reader],plan['closure_plan'])
    print(json.dumps(dict(status='queued_complete_original_grid_recovery_diagnostics_and_readback',
        launches=[producer,reader,closer],gpu=False,native_sampling_restarted=False)),flush=True)


if __name__=='__main__':main()
