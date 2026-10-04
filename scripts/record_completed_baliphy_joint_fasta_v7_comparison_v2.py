#!/usr/bin/env python3
"""Close completed native output while retaining the original reporting wrapper failure."""
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import live_group,sha
from baliphy_joint_node_logger_v3 import validate_frame
from independent_native_ancestral_alignment import fasta_records,tree_labels
from read_baliphy_scalar_json_v6b import compare_tsv
from reference_measurement_union_sources import bind,verify


def main():
    output=Path('metadata/baliphy_joint_fasta_v7_full_comparison_review_20261004_v2.json')
    assert not output.exists()
    pp=Path('metadata/baliphy_joint_fasta_v7_full_comparison_plan_20261004_v1.json')
    ep=Path('metadata/baliphy_joint_fasta_v7_full_comparison_execution_20261004_v1.json')
    tp=Path('metadata/baliphy_joint_fasta_v7_full_comparison_original_tool_payloads_20261004_v1.json')
    plan,execution,payload=[json.loads(p.read_text()) for p in [pp,ep,tp]]
    assert payload['original_tool_session_id']==payload['initial']['session_id']==58381
    assert payload['terminal']['exit_code']==1
    assert execution['status']=='failed_stage_retained' and execution['exit_code']==1
    assert not execution['timed_out'] and execution['receipt_sha256'] is None
    inv=execution['invocation_id'];assert inv in payload['initial']['output']
    unit='fungal-joint-fasta-v7-full-comparison-20261004-v1.service'
    records=[json.loads(line) for line in subprocess.check_output(
        ['journalctl','--user','-u',unit,'-o','json','--no-pager'],text=True).splitlines()]
    records=[r for r in records if inv in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
    exact=[r for r in records if r.get('_PID')==str(execution['wrapper']['pid'])
           and r.get('_CMDLINE')==' '.join(execution['wrapper']['cmdline'])]
    assert len(exact)==2
    assert json.loads(exact[0]['MESSAGE'])==dict(original_wrapper=execution['wrapper'],invocation_id=inv)
    terminal={k:v for k,v in execution.items() if k not in ['source_hashes','artifacts','command','wrapper','child','scope']}
    assert json.loads(exact[1]['MESSAGE'])==terminal
    starts=[r for r in records if r.get('USER_INVOCATION_ID')==inv and 'Started ' in r.get('MESSAGE','')]
    ends=[r for r in records if r.get('USER_INVOCATION_ID')==inv and r.get('CPU_USAGE_NSEC')]
    assert len(starts)==len(ends)==1
    journal=ep.with_suffix('')/'original-invocation-journal.jsonl'
    with journal.open('x') as f:
        for r in records:f.write(json.dumps(r,sort_keys=True)+'\n')
    pins=dict(plan['pins'])
    for mapping in [execution['source_hashes'],execution['artifacts']]:
        for path,h in mapping.items():bind(pins,path,h)
    verify(pins)
    stderr=Path(ep.with_suffix(''),'stderr.log').read_text()
    assert 'FileNotFoundError' in stderr and 'C1.log.column-map.json' in stderr
    root=Path(plan['output'])/'native/attempt-0001'
    native_path=root/'receipt.json';native=json.loads(native_path.read_text())
    identity=json.loads((root/'process.json').read_text())
    assert native['exit_code']==0 and native['status']=='exited_zero_pending_scientific_validation'
    assert not live_group(identity['pgid'])
    assert json.loads((root.parent/'configuration.json').read_text())==plan['native_config']
    assert json.loads((root/'command.json').read_text())==identity['command']==plan['native_config']['command']
    for name,h in native['artifacts'].items():bind(pins,root/name,h)
    directory=root/'independent-chain-1'
    comparison=compare_tsv(directory)
    assert comparison['rows']==21 and comparison['mapped_values_compared']==903
    labels,tips=tree_labels((directory/'runtime-tree.nwk').read_text())
    assert len(labels)==1243 and len(tips)==622
    records=[json.loads(line) for line in (directory/'C1.P1.site-property-samples.jsonl').read_text().splitlines()]
    assert [r['iter'] for r in records]==[0,10,20]
    frames=[]
    for record in records:
        sequences=fasta_records(record['alignmentLines']);assert set(sequences)==labels
        frames.append(validate_frame(record,record['iter'],sequences,tips))
    before=Path('results/ancestral/scalar-native-stack-comparison-20261004-v2/independent-chain-1')
    prefix={}
    for name in ['C1.log','C1.log.json','C1.log.column-map.json','runtime-tree.nwk',
                 'C1.P1.fastas','C1.P1.site-property-samples.jsonl']:
        old=before/name;new=directory/name;a=old.read_bytes();b=new.read_bytes()
        assert a and b.startswith(a),(name,'nonempty original64MiB prefix differs')
        prefix[name]=dict(original_bytes=len(a),new_bytes=len(b),nonempty_prefix_identical=True)
        bind(pins,old)
    for path in [Path(__file__),pp,ep,tp,journal,native_path,root.parent/'configuration.json']:
        bind(pins,path)
    verify(pins)
    result=dict(status='completed_candidate_full_input_horizon_pending_full_failure_grid_qualification',
        checked_utc=datetime.now(timezone.utc).isoformat(),plan_sha256=sha(pp),
        native_exit_code=0,native_status=native['status'],native_receipt=str(native_path),
        native_elapsed_seconds=native['elapsed_seconds'],planned_horizon_completed=True,
        native_caps=plan['native_caps'],scalar_readback=comparison,joint_frames=frames,
        nonempty_original_64_mib_prefix_comparisons=prefix,
        original_wrapper_reporting_failure=dict(original_tool_session_id=58381,
            original_tool_exit_code=1,invocation_id=inv,wrapper=execution['wrapper'],
            whole_wrapper_initial_and_terminal_payloads_matched=True,manager_start_records=1,
            manager_completion_records=1,error_type='FileNotFoundError',
            missing_original8MiB_diagnostic_file='C1.log.column-map.json',
            original_v1_source_and_output_namespace_preserved=True),
        child_peak_rss_bytes_as_recorded_by_wrapper=execution['child_peak_rss_bytes'],
        cgroup_manager_memory_is_not_native_peak=True,
        source_hashes=pins,full_24_original_failure_grid_qualified=False,validated_production_repair=False,
        original_jobs_restarted=False,installed_software_changed=False,
        scientific_eligibility=False,posterior_qualified=False,gpu=False,
        scope='Original standalone native exits0 with21scalar rows and3fulljointframes at original '
              '8MiB stack and other native caps. Wrapper originally exits1 on missing older8MiB '
              'diagnostic column-map during reporting; that exact failure/source/tool/journal is '
              'preserved. Fresh readonlyV2 uses nonempty original64MiB evidence, verifies all '
              'native artifacts,903scalar values and3primary full-node frames, without rerunning '
              'any native attempt. Independent full source/tree/projection replay is still separate; '
              'all24failure roles and adequate posterior remain unqualified.')
    with output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','joint_frames','nonempty_original_64_mib_prefix_comparisons']},indent=2))


if __name__=='__main__':main()
