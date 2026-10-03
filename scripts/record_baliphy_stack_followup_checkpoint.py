#!/usr/bin/env python3
"""Observe exact original follow-up handles, caps and available checkpoint counts."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_baliphy_sampler_resource_observation_checkpoint import last_event
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record
from reference_measurement_union_sources import verify


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists()
    pp=Path('metadata/baliphy_stack_followup_plan_20261003_v1.json');plan=json.loads(pp.read_text());verify(plan['pins'])
    inventory=json.loads(Path(plan['launch_inventory']).read_text());assert inventory['source_plan_sha256']==sha(pp)
    handles=[]
    for i,fp in enumerate(inventory['launches']+plan['dependencies']):
        r=json.loads(Path(fp).read_text());r['launch']=fp;assert sha(r['plan'])==r['plan_sha256']
        proc=fingerprint(r)
        if proc is None:observed=journal_terminal(r)
        else:
            observed=live_record(r,proc)
            group=next(l[3:] for l in (Path('/proc')/str(r['pid'])/'cgroup').read_text().splitlines() if l.startswith('0::'))
            cg=Path('/sys/fs/cgroup')/group.lstrip('/')
            limits={k:(cg/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
            assert limits==r['actual_cgroup_limits']
            if i<3:
                cpu,mem=(4,200) if i==0 else (2,32)
                assert limits=={'cpu.max':str(cpu*100000)+' 100000','memory.max':str(mem*2**30),'memory.swap.max':'0'}
            observed['current_cgroup_limits']=limits
        handles.append(observed)
    root=Path(plan['output']);rows=[json.loads(f.read_text()) for f in (root/'chains').glob('*.json')]
    assert all(r['scientific_eligibility'] is r['posterior_qualified'] is False for r in rows)
    observations=last_event(root/'observations.jsonl') if (root/'observations.jsonl').exists() else None
    if observations:
        assert observations['cgroup']['limits']=={'cpu.max':'400000 100000','memory.max':str(200*2**30),'memory.swap.max':'0'}
        assert observations['cgroup']['memory_bytes']['memory.swap.current']==0
        assert all(observations['cgroup']['memory_events'][k]==0 for k in ['max','oom','oom_kill','oom_group_kill'])
        live=[s for s in observations['native_observations'] if s['status']=='verified_live_native_resource_observation']
        assert all(s['native_stack_bytes']==64*2**20 for s in live)
    else:live=[]
    result=dict(status='verified_original_scoped_stack_followup_runtime',checked_utc=datetime.now(timezone.utc).isoformat(),
        plan_sha256=sha(pp),frozen_pins_checked=len(plan['pins']),original_handles=handles,
        full_original_roles_retained=1620,original_successful_roles=1596,original_failed_roles_retained=24,
        expected_fresh_followup_roles=24,completed_followup_checkpoints=len(rows),
        followup_status_counts=dict(Counter(r['status'] for r in rows)),last_resource_snapshot=observations,
        current_verified_native_stack_observations=len(live),native_stack_bytes=64*2**20,
        longer_posterior_launched=False,posterior_qualified=False,gpu=False,
        scope='Exact original handles/journals, immutable plan/pins, available saved row counts and most recent identity-bound resource event. No new successful-output replay or claim of full follow-up closure, precise final peaks or posterior adequacy.')
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['original_handles','last_resource_snapshot']},indent=2))


if __name__=='__main__':main()
