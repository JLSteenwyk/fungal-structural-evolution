#!/usr/bin/env python3
"""Independently reread the full candidate with the existing native/topology/scalar/joint readers."""
import argparse
import copy
from datetime import datetime,timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_joint_fasta_lines_v7 import reverse
from baliphy_joint_sampler_scalar_v6 import inspect,SUCCESS
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['plan','source-receipt','transport','output','receipt']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists()
    plan=json.loads(a.plan.read_text());source=json.loads(a.source_receipt.read_text())
    transport=json.loads(a.transport.read_text())
    assert source['status']=='completed_candidate_full_input_horizon_pending_full_failure_grid_qualification'
    assert source['native_exit_code']==0 and source['planned_horizon_completed']
    assert source['plan_sha256']==sha(a.plan)
    assert transport['validation_sha256']==sha(a.source_receipt)
    assert transport['original_tool_terminal_exit_code']==0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    pins=dict(transport['source_hashes']);verify(pins)
    original_plan_path=Path('metadata/baliphy_scalar_v6_sampler_execution_plan_20261004_v1.json')
    original_plan=json.loads(original_plan_path.read_text())
    jobs_path=Path(original_plan['jobs']);jobs=json.loads(jobs_path.read_text());assert len(jobs)==1620
    config=plan['native_config'];command=config['command']
    seed=int(command[command.index('--seed')+1]);original_program=plan['original_model']
    matches=[j for j in jobs if j['chain']['seed']==seed and j['chain']['program']==original_program]
    assert len(matches)==1
    original=matches[0];job=copy.deepcopy(original)
    model=Path(plan['diagnostic_model']);assert reverse(model.read_text())==Path(original_program).read_text()
    job['chain'].update(chain_id=original['chain']['chain_id']+'-joint-fasta-v7-full-comparison',
                        program=str(model),program_sha256=sha(model))
    job.update(config=config,source_seed=original['chain']['seed'],source_chain_id=original['chain']['chain_id'])
    root=a.output.resolve();root.mkdir(exist_ok=False)
    native_receipt=Path(source['native_receipt'])
    row=inspect(job,native_receipt,sha(a.plan),original_plan['mapping'],root/'frames')
    assert row['status']==SUCCESS and row['scalar_integrity_accepted']
    assert row['saved_alignments']==3 and row['candidate_frames']==12
    assert row['scalar_v6_audit']['mapped_values_compared']==903
    assert len(row['joint_frames'])==3 and all(f['native_tips']==622 for f in row['joint_frames'])
    frames=row['joint_frames']
    assert all(f['native_nodes']==1243 and f['native_ancestors']==621 for f in frames)
    directory=native_receipt.parent/'independent-chain-1'
    previous=Path('results/ancestral/scalar-native-stack-comparison-20261004-v2/independent-chain-1')
    prefix={}
    for name in ['C1.log','C1.log.json','C1.log.column-map.json','runtime-tree.nwk',
                 'C1.P1.fastas','C1.P1.site-property-samples.jsonl']:
        old=previous/name;new=directory/name;before=old.read_bytes();after=new.read_bytes()
        assert before and after.startswith(before),(name,'nonempty original64MiB prefix differs')
        prefix[name]=dict(original_bytes=len(before),new_bytes=len(after),nonempty_prefix_identical=True)
        bind(pins,old)
    export=root/'disposition.json'
    with export.open('x') as f:json.dump(row,f,indent=2,allow_nan=False);f.write('\n')
    for path in [Path(__file__),a.plan,a.source_receipt,a.transport,original_plan_path,jobs_path,
                 Path(original_plan['mapping']),model,export,
                 Path('scripts/baliphy_joint_sampler_scalar_v6.py'),
                 Path('scripts/baliphy_joint_sampler_qualification_v3.py'),
                 Path('scripts/independent_joint_ancestral_frames.py')]:bind(pins,path)
    for frame in frames:bind(pins,frame['projection_array'],frame['projection_array_sha256'])
    verify(pins)
    result=dict(status='passed_independent_full_v7_candidate_comparison_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(),source_receipt_sha256=sha(a.source_receipt),
        original_source_chain_id=original['chain']['chain_id'],native_receipt=str(native_receipt),
        scalar_rows=21,mapped_scalar_values=903,joint_frames=3,candidate_frames=12,
        native_tips=622,native_ancestors=621,
        ancestral_residue_category_pairs=sum(f['ancestral_pairs'] for f in frames),
        tip_residue_category_pairs=sum(f['tip_pairs'] for f in frames),
        nonempty_original_64_mib_prefix_comparisons=prefix,
        original_8_mib_empty_files_excluded_from_compatibility_evidence=True,
        disposition=str(export),source_hashes=pins,
        full_24_failure_grid_qualified=False,scientific_eligibility=False,posterior_qualified=False,
        original_jobs_restarted=False,installed_software_changed=False,gpu=False,
        scope='Existing independent native/tree/branch/FASTA/scalar/category readers check all '
              'three full1243-node frames and mapped source candidates. Original native command '
              'and artifact custody are reconstructed, and all ten saved64MiB scalar rows plus '
              'the nonempty first joint/FASTA frame match exactly. Only output integrity for '
              'one full failed input is qualified, not all24original failures or posterior mixing.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','nonempty_original_64_mib_prefix_comparisons']},indent=2))


if __name__=='__main__':main()
