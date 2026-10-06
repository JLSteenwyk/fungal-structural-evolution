#!/usr/bin/env python3
"""Launch the v2 full comparison export from an already verified closure."""
import argparse
import json
from pathlib import Path
import shutil
import sys

from full_whole_protein_comparison_sources import SUMMARY_FIELDS
from launch_full_triad_sequence_geometry_followup import create,launch
from screen_duplication_alignment_reuse import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text())
    for path,digest in plan['pins'].items():assert sha(path)==digest,path
    closure=json.loads(Path(plan['optimization_completion']).read_text())
    assert closure['status']=='complete_verified_full_whole_protein_optimization'
    assert closure['scientific_eligibility'] is False
    assert closure['exact_process_journals_checked']==plan['expected_process_journals']
    assert closure['historical_source_hashes']==len(plan['historical_source_hashes'])
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    fixture=json.loads(Path(plan['fixture_validation']).read_text())
    assert fixture['status']=='passed_full_whole_protein_comparison_software_contracts'
    assert fixture['rejected_false_exports']==26
    assert fixture['checked_reused_tree_checkpoints']==1
    assert fixture['completed_restart_refused'] is True
    root=Path(plan['output']);assert not root.exists()
    prefix='metadata/full_whole_protein_comparisons_v2'
    producer=launch('full-whole-protein-comparisons-v2',prefix+'_wait_plan_20261006.json',
        [sys.executable,'scripts/export_full_whole_protein_comparisons.py','--plan',str(args.plan)],[],str(args.plan))
    reader=launch('full-whole-protein-comparisons-v2-readback',prefix+'_readback_wait_plan_20261006.json',
        [sys.executable,'scripts/readback_full_whole_protein_comparisons.py','--plan',str(args.plan),'--output',str(root/'readback.json')],[producer],str(args.plan))
    completion=prefix+'_completion_plan_20261006.json'
    create(completion,dict(source_plan=str(args.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status='complete_full_whole_protein_comparisons_pending_independent_readback',reader_status='passed_full_whole_protein_comparisons_source_selection_arithmetic_readback',
        completed_status='complete_verified_full_whole_protein_comparisons',summary_fields=SUMMARY_FIELDS,launches=[producer,reader],
        pins={str(path):sha(path) for path in [args.plan,producer,reader,'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=prefix+'_completed_20261006.json',scope=plan['scope']))
    closer=launch('full-whole-protein-comparisons-v2-closure',prefix+'_closure_wait_plan_20261006.json',
        [sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',completion],[producer,reader],completion)
    create(prefix+'_launches_20261006.json',dict(status='queued_full_whole_protein_comparison_integration_from_verified_v2_closure',
        source_plan=str(args.plan),source_plan_sha256=sha(args.plan),launches=[producer,reader,closer],scope=plan['scope']))


if __name__=='__main__':main()
