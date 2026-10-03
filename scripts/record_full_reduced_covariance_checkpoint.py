#!/usr/bin/env python3
"""Verify original retained-basis queue handles, source gates and live caps."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record
from reference_measurement_union_sources import verify


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    path=Path('metadata/full_reduced_covariance_qualification_plan_20261003_v1.json');plan=json.loads(path.read_text())
    verify(plan['pins']);inventory=json.loads(Path(plan['launch_inventory']).read_text())
    assert inventory['source_plan_sha256']==sha(path)
    handles=[]
    for index,launch in enumerate(inventory['launches']+plan['dependencies']):
        row=json.loads(Path(launch).read_text());row['launch']=launch
        assert sha(row['plan'])==row['plan_sha256'];process=fingerprint(row)
        if process is None:observed=journal_terminal(row)
        else:
            observed=live_record(row,process)
            group=next(s[3:] for s in Path(f'/proc/{process.pid}/cgroup').read_text().splitlines() if s.startswith('0::'))
            root=Path('/sys/fs/cgroup')/group.lstrip('/')
            limits={k:(root/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
            memory=32 if index==3 else 16
            assert limits==row['actual_cgroup_limits']=={'cpu.max':'200000 100000','memory.max':str(memory*2**30),'memory.swap.max':'0'}
            observed['current_cgroup_limits']=limits
        handles.append(observed)
    gates={}
    for key in ['qualification_completion','exact_completion','completion']:
        fp=Path(plan[key]);entry=dict(path=str(fp),closure_present=fp.exists())
        if fp.exists():
            value=json.loads(fp.read_text());archive=Path(value['full_hash_archive'])
            assert sha(archive)==value['full_hash_archive_sha256']
            entry.update(status=value['status'],sha256=sha(fp),full_hash_archive_sha256=value['full_hash_archive_sha256'],
                declared_source_bindings=value['bound_source_hashes'],original_journals=value['exact_process_journals_checked'])
        gates[key]=entry
    result=dict(status='verified_original_complete_retained_covariance_qualification_queue',checked_utc=datetime.now(timezone.utc).isoformat(),
        plan_sha256=sha(path),frozen_pins_checked=len(plan['pins']),original_handles=handles,completion_gates=gates,
        expected=plan['expected'],production_fitting_launched=False,gpu=False,all_eight_aims_incomplete=True,
        scope='Five original PID/create/command or invocation-linked terminal journals and live caps; frozen scripts/plans/software proofs and compact archive hashes checked. Broader inherited files not rehashed here. Missing original arithmetic closure keeps new native work queued; no restart, fit or biological acceptance.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='original_handles'}))


if __name__=='__main__':main()
