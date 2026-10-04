#!/usr/bin/env python3
"""Capture a closed prefix from a live candidate; compare nonempty original64MiB output."""
from datetime import datetime,timezone
import json
from pathlib import Path

import psutil

from ancestral_chain_attempt import sha
from baliphy_joint_node_logger_v3 import validate_frame
from independent_native_ancestral_alignment import fasta_records,tree_labels
from read_baliphy_scalar_json_v6b import compare_tsv
from reference_measurement_union_sources import bind,verify


def complete_lines(path):
    data=path.read_bytes()
    return data[:data.rfind(b'\n')+1].splitlines(keepends=True)


def main():
    output=Path('metadata/baliphy_joint_fasta_v7_prefix_checkpoint_20261004_v1.json')
    assert not output.exists()
    pp=Path('metadata/baliphy_joint_fasta_v7_full_comparison_plan_20261004_v1.json')
    plan=json.loads(pp.read_text());verify(plan['pins'])
    pins={}
    base=Path(plan['output'])/'native/attempt-0001'
    identity_path=base/'process.json';identity=json.loads(identity_path.read_text())
    native=psutil.Process(identity['pid'])
    assert native.create_time()==identity['created']
    command=identity['command'];command=command[command.index('--')+1:]
    assert native.cmdline()==command
    caps=plan['native_caps']
    assert native.rlimit(psutil.RLIMIT_AS)==(caps['address_space_bytes'],)*2
    assert native.rlimit(psutil.RLIMIT_CPU)==(caps['cpu_seconds'],)*2
    assert native.rlimit(psutil.RLIMIT_FSIZE)==(caps['per_file_bytes'],)*2
    assert native.rlimit(psutil.RLIMIT_STACK)==(8*2**20,-1)
    before_path=Path('metadata/scalar_native_stack_comparison_terminal_review_20261004_v3.json')
    before_review=json.loads(before_path.read_text());verify(before_review['source_hashes'])
    assert before_review['native_signal_from_original_debugger_text']=='SIGKILL'
    before=Path('results/ancestral/scalar-native-stack-comparison-20261004-v2/independent-chain-1')
    current=base/'independent-chain-1'
    snapshot=Path('data/software_audits/baliphy-joint-fasta-v7-prefix-observation-20261004-v1')
    snapshot.mkdir(exist_ok=False)
    tsv=complete_lines(current/'C1.log');scalar=complete_lines(current/'C1.log.json')
    row_count=min(len(tsv)-1,len(scalar)-1)
    assert 1<=row_count<=10
    for name,rows in [('C1.log',tsv),('C1.log.json',scalar)]:
        captured=b''.join(rows[:row_count+1])
        assert len(captured)>0 and (before/name).read_bytes().startswith(captured)
        (snapshot/name).write_bytes(captured);bind(pins,before/name)
    for name in ['C1.log.column-map.json','runtime-tree.nwk']:
        captured=(current/name).read_bytes();assert captured==(before/name).read_bytes()
        (snapshot/name).write_bytes(captured);bind(pins,before/name)
    joint_lines=complete_lines(current/'C1.P1.site-property-samples.jsonl')
    assert joint_lines
    captured=joint_lines[0]
    assert captured==complete_lines(before/'C1.P1.site-property-samples.jsonl')[0]
    (snapshot/'C1.P1.site-property-samples.jsonl').write_bytes(captured)
    bind(pins,before/'C1.P1.site-property-samples.jsonl')
    scalar_comparison=compare_tsv(snapshot)
    assert scalar_comparison['rows']==row_count
    frame=json.loads(captured);assert frame['iter']==0
    labels,tips=tree_labels((snapshot/'runtime-tree.nwk').read_text())
    assert len(labels)==1243 and len(tips)==622
    sequences=fasta_records(frame['alignmentLines']);assert set(sequences)==labels
    joint=validate_frame(frame,0,sequences,tips)
    original8=Path('results/ancestral/scalar-native-signal-debugger-20261004-v1/independent-chain-1')
    assert (original8/'C1.log').stat().st_size==(original8/'C1.P1.site-property-samples.jsonl').stat().st_size==0
    for p in [Path(__file__),pp,identity_path,before_path,*snapshot.iterdir(),
              original8/'C1.log',original8/'C1.P1.site-property-samples.jsonl']:bind(pins,p)
    verify(pins)
    assert native.create_time()==identity['created'] and native.cmdline()==command
    result=dict(status='verified_live_original_limit_candidate_with_closed_identical_saved_prefix',
        checked_utc=datetime.now(timezone.utc).isoformat(),original_tool_session_id=58381,
        native=dict(pid=native.pid,created=identity['created'],cmdline=command,
                    status=native.status(),cpu_seconds=native.cpu_times().user+native.cpu_times().system,
                    address_space=native.rlimit(psutil.RLIMIT_AS),cpu_limit=native.rlimit(psutil.RLIMIT_CPU),
                    file_limit=native.rlimit(psutil.RLIMIT_FSIZE),stack_limit=native.rlimit(psutil.RLIMIT_STACK)),
        snapshot_root=str(snapshot),scalar_readback=scalar_comparison,joint_frame_readback=joint,
        captured_first_joint_frame_bytes=len(captured),
        nonempty_scalar_and_joint_prefix_byte_identical_to_original_64_mib_run=True,
        original_8_mib_uninstrumented_saved_scalar_rows=0,
        original_8_mib_uninstrumented_saved_joint_frames=0,
        earlier_trace_review_v1_iteration_comparison_needs_this_qualification=True,
        full_horizon_verified=False,all_24_failures_repaired=False,
        source_hashes=pins,original_jobs_restarted=False,scientific_eligibility=False,posterior_qualified=False,
        scope='Live original8MiB stack candidate produces the first full joint frame and scalar rows. '
              'A closed captured prefix matches actual nonempty output of the original64MiB run and '
              'is independently reread. Original8MiB uninstrumented debugger files are empty, '
              'so their vacuous prefix comparisons provide no compatibility evidence. The phase '
              'diagnostic also failed before its first joint frame; its V1 statement about later '
              'uninstrumented iterations applies to the separate64MiB run, not the original8MiB run. '
              'No completed20-iteration horizon, all24-failure repair or posterior acceptance is inferred.')
    with output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','native','joint_frame_readback']},indent=2))


if __name__=='__main__':main()
