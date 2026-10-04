#!/usr/bin/env python3
"""Run one fresh candidate attempt and independently check a completed output horizon."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

from ancestral_chain_attempt import run_attempt,sha
from baliphy_joint_node_logger_v3 import validate_frame
from independent_native_ancestral_alignment import fasta_records,tree_labels
from read_baliphy_scalar_json_v6b import compare_tsv
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists()
    plan=json.loads(a.plan.read_text());pins=dict(plan['pins']);verify(pins)
    root=Path(plan['output']).resolve();root.mkdir(exist_ok=False)
    receipt_path=run_attempt(root/'native',plan['native_config'])
    native=json.loads(receipt_path.read_text());comparison=None;frames=[];prefix={}
    if native['exit_code']==0:
        directory=receipt_path.parent/'independent-chain-1'
        comparison=compare_tsv(directory)
        assert comparison['rows']==21 and comparison['mapped_values_compared']==903
        labels,tips=tree_labels((directory/'runtime-tree.nwk').read_text())
        assert len(tips)==622 and len(labels)==1243
        records=[json.loads(line) for line in (directory/'C1.P1.site-property-samples.jsonl').read_text().splitlines()]
        assert [r['iter'] for r in records]==[0,10,20]
        for record in records:
            sequences=fasta_records(record['alignmentLines']);assert set(sequences)==labels
            frames.append(validate_frame(record,record['iter'],sequences,tips))
        before=Path(plan['original_debugger_directory'])
        for name in ['C1.log','C1.log.json','C1.log.column-map.json','runtime-tree.nwk',
                     'C1.P1.fastas','C1.P1.site-property-samples.jsonl']:
            old=before/name;new=directory/name
            bind(pins,old)
            prefix[name]=dict(original_bytes=old.stat().st_size,new_bytes=new.stat().st_size,
                              original_complete_bytes_are_prefix=new.read_bytes().startswith(old.read_bytes()))
    for path in [a.plan,Path(__file__),receipt_path,root/'native/configuration.json']:bind(pins,path)
    pins.update({str(receipt_path.parent/name):digest for name,digest in native['artifacts'].items()})
    verify(pins)
    result=dict(status='completed_candidate_full_input_horizon_pending_full_failure_grid_qualification'
        if native['exit_code']==0 else 'failed_separate_candidate_full_input_comparison',
        checked_utc=datetime.now(timezone.utc).isoformat(),plan_sha256=sha(a.plan),
        native_exit_code=native['exit_code'],native_status=native['status'],native_receipt=str(receipt_path),
        native_elapsed_seconds=native['elapsed_seconds'],planned_horizon_completed=native['exit_code']==0,
        native_caps=plan['native_caps'],scalar_readback=comparison,joint_frames=frames,
        original_saved_output_prefix_comparisons=prefix,source_hashes=pins,
        full_24_original_failure_grid_qualified=False,validated_production_repair=False,
        original_jobs_restarted=False,installed_software_changed=False,
        scientific_eligibility=False,posterior_qualified=False,gpu=False,scope=plan['scope'])
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','joint_frames']},indent=2))


if __name__=='__main__':main()
