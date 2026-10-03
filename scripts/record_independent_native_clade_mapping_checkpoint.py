#!/usr/bin/env python3
"""Observe exact native replay handles, resource caps and unclosed progress."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,default=Path('metadata/independent_native_clade_mapping_plan_20261002.json'))
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text())
    pins={str(Path(__file__)):sha(__file__),str(a.plan):sha(a.plan),**plan['pins']}
    assert all(sha(k)==v for k,v in pins.items())
    inventory=json.loads(Path(plan['launch_inventory']).read_text());assert inventory['source_plan_sha256']==sha(a.plan)
    handles=[]
    for path in inventory['launches']:
        record=json.loads(Path(path).read_text());record['launch']=path
        assert sha(record['plan'])==record['plan_sha256'];proc=fingerprint(record)
        if proc is None:observation=journal_terminal(record)
        else:
            observation=live_record(record,proc)
            try:
                lines=(Path('/proc')/str(record['pid'])/'cgroup').read_text().splitlines()
                unified=[line[3:] for line in lines if line.startswith('0::')];assert len(unified)==1
                root=Path('/sys/fs/cgroup')/unified[0].lstrip('/')
                limits={k:(root/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
                assert limits==record['actual_cgroup_limits']=={'cpu.max':'200000 100000','memory.max':str(32*2**30),'memory.swap.max':'0'}
                observation['current_cgroup_limits']=limits
                observation['observed_memory_bytes']={k:(root/k).read_text().strip() for k in ['memory.current','memory.peak','memory.swap.current']}
            except FileNotFoundError:
                assert fingerprint(record) is None;observation=journal_terminal(record)
        handles.append(observation)
    root=Path(plan['output']);checked=failed=node_pairs=branches=candidates=frames=roots=0
    checkpoints=sorted((root/'chains').glob('*.json'))
    for path in checkpoints:
        row=json.loads(path.read_text());assert row['chain_id']==path.stem and row['stage']['plan_sha256']==sha(a.plan)
        result=row['result'];assert result['scientific_eligibility'] is False
        if result['status']=='unresolved_failed_native_chain_retained':
            assert result['node_rows']==result['candidate_rows']==[];failed+=1
        else:
            assert result['status']=='all_source_runtime_rooted_clades_and_candidate_mappings_checked_conditional_on_root'
            assert result['biological_root_accepted'] is False;checked+=1
            node_pairs+=len(result['node_rows']);branches+=sum(not r['is_root'] for r in result['node_rows'])
            candidates+=len(result['candidate_rows']);frames+=result['candidate_frame_mappings']
            roots+=sum(r['assumed_root'] for r in result['candidate_rows'])
    completion=None;cp=Path(plan['completion'])
    if cp.exists():
        completion=json.loads(cp.read_text());assert completion['status']=='complete_verified_full_independent_native_clade_mapping'
        assert sha(completion['full_hash_archive'])==completion['full_hash_archive_sha256']
    result=dict(status='verified_original_native_clade_mapping_runtime_and_unclosed_checkpoint_observation',
        checked_utc=datetime.now(timezone.utc).isoformat(),source_hashes=pins,original_handles=handles,expected=plan['expected'],
        unclosed_producer_checkpoints=dict(chains=len(checkpoints),checked_chains=checked,failed_chains=failed,
            node_pairs=node_pairs,nonroot_branch_pairs=branches,candidate_node_pairs=candidates,candidate_frame_mappings=frames,assumed_root_candidates=roots),
        accounting_closure=completion,gpu_inference_paused=True,all_eight_aims_incomplete=True,
        scope='Exact original PID/create/command or actual completion journal, small frozen pins and live caps observed. '
              'Checkpoint totals are unclosed progress, not full serialized/source/journal closure or posterior acceptance. '
              'Compact closure archive hash, if present, is checked; this observer does not repeat full topology reconstruction or large-data hashing.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','original_handles','accounting_closure']}))


if __name__=='__main__':main()
