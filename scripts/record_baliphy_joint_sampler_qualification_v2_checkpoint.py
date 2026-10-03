#!/usr/bin/env python3
"""Observe original joint-grid controllers and all four dependency closures."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_project_runtime_checkpoint_v4 import fingerprint,journal_terminal,live_record
from run_baliphy_reference_preflight import verify


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    path=Path('metadata/baliphy_joint_sampler_qualification_v2_plan_20261003.json');plan=json.loads(path.read_text());verify(plan)
    inventory=json.loads(Path(plan['launch_inventory']).read_text());assert inventory['source_plan_sha256']==sha(path)
    handles=[]
    for index,launch in enumerate(inventory['launches']+plan['dependencies']):
        record=json.loads(Path(launch).read_text());record['launch']=launch
        assert sha(record['plan'])==record['plan_sha256'];process=fingerprint(record)
        if process is None:observed=journal_terminal(record)
        else:
            observed=live_record(record,process)
            group=next(line[3:] for line in (Path('/proc')/str(record['pid'])/'cgroup').read_text().splitlines() if line.startswith('0::'))
            root=Path('/sys/fs/cgroup')/group.lstrip('/')
            limits={k:(root/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
            assert limits==record['actual_cgroup_limits']
            cpus,memory=[(16,200),(2,32),(2,32),(2,32),(2,32),(2,4),(2,16)][index]
            assert limits=={'cpu.max':str(cpus*100000)+' 100000','memory.max':str(memory*2**30),'memory.swap.max':'0'}
            observed['current_cgroup_limits']=limits
        handles.append(observed)
    rows=[json.loads(p.read_text()) for p in (Path(plan['output'])/'chains').glob('*.json')]
    assert all(r['scientific_eligibility'] is r['posterior_qualified'] is False for r in rows)
    result=dict(status='verified_original_full_joint_sampler_v2_runtime',checked_utc=datetime.now(timezone.utc).isoformat(),
        plan_sha256=sha(path),frozen_pins_checked=len(plan['pins']),original_handles=handles,
        unclosed_joint_role_checkpoints=len(rows),joint_status_counts=dict(Counter(r['status'] for r in rows)),
        expected_native_roles=1620,iterations_per_role=20,longer_posterior_horizon_launched=False,
        posterior_qualified=False,gpu=False,scope='Exact seven original handles/journals and live limits, all immutable pins and unclosed row counts; not native output replay, full closure, longer resource or posterior qualification.')
    with a.output.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='original_handles'}))


if __name__=='__main__':main()
