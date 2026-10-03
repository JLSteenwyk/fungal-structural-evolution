#!/usr/bin/env python3
"""Live native procfs checks and full-grid missing-observation accounting."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import time
from unittest.mock import Mock, patch

import psutil

from ancestral_chain_attempt import run_attempt, sha, write_json
from baliphy_sampler_resource_observation import describe, probe_native, parse_status, parse_limits, finish_roles, summarize
from run_baliphy_reference_preflight import verify


def rejected(action):
    try:action()
    except (AssertionError,KeyError,ValueError):return
    raise AssertionError('Unsafe or malformed memory observation accepted')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--receipt',type=Path,required=True);args=parser.parse_args()
    root=args.output.resolve();root.mkdir(exist_ok=False)
    native_path=Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json')
    native=json.loads(native_path.read_text());verify(native)
    jobs=json.loads(Path(native['jobs']).read_text());assert len(jobs)==1620
    status='VmPeak:\t2048 kB\nVmSize:\t1024 kB\nVmHWM:\t768 kB\nVmRSS:\t512 kB\nVmSwap:\t0 kB\n'
    expected={name:value*1024 for name,value in [('VmPeak',2048),('VmSize',1024),('VmHWM',768),('VmRSS',512),('VmSwap',0)]}
    assert parse_status(status)==expected;negatives=[]
    for name,text in [('duplicate_memory_field',status+'VmRSS: 1 kB\n'),('missing_memory_field',status.replace('VmPeak:\t2048 kB\n','')),
                      ('negative_memory',status.replace('512 kB','-1 kB')),('wrong_memory_unit',status.replace('512 kB','512 MB'))]:
        rejected(lambda:parse_status(text));negatives.append(name)
    job=jobs[0];cmd=job['config']['command'];fake=root/'proc'/'123';fake.mkdir(parents=True)
    (fake/'status').write_text(status);(fake/'cgroup').write_text('0::/software-fixture\n')
    limits='Max address space  %d  %d bytes\nMax cpu time %d %d seconds\nMax file size %d %d bytes\n' % (
        job['memory_reservation_bytes'],job['memory_reservation_bytes'],job['config']['timeout_seconds'],
        job['config']['timeout_seconds'],2*2**30,2*2**30)
    (fake/'limits').write_text(limits);assert parse_limits(limits,job)['Max address space']==job['memory_reservation_bytes']
    descriptor=dict(chain_id=job['chain']['chain_id'],pid=123,created=100.0,identity_sha256='software-only')
    process=Mock();process.create_time.return_value=100.0;process.status.return_value=psutil.STATUS_RUNNING
    process.cmdline.return_value=cmd[cmd.index('--')+1:];process.cpu_times.return_value=(1.,2.)
    with patch('baliphy_sampler_resource_observation.psutil.Process',return_value=process):
        observed=probe_native(job,descriptor,'/software-fixture',root/'proc')
        assert observed['reported_memory_bytes']==expected and observed['cpu_seconds']==3
        process.create_time.return_value=101.0
        assert probe_native(job,descriptor,'/software-fixture',root/'proc')['status']=='original_process_unavailable_at_observation'
        process.create_time.return_value=100.0;process.cmdline.return_value=['unrelated-worker']
        rejected(lambda:probe_native(job,descriptor,'/software-fixture',root/'proc'));negatives.append('changed_native_command')
        process.cmdline.return_value=[]
        assert probe_native(job,descriptor,'/software-fixture',root/'proc')['status']=='native_command_unavailable_or_transitioning_at_observation'
        process.cmdline.return_value=cmd
        assert probe_native(job,descriptor,'/software-fixture',root/'proc')['status']=='limit_wrapper_before_native_exec_at_observation'
        process.cmdline.return_value=cmd[cmd.index('--')+1:]
        rejected(lambda:probe_native(job,descriptor,'/different-cgroup',root/'proc'));negatives.append('different_native_cgroup')
        (fake/'limits').write_text(limits.replace(str(job['memory_reservation_bytes']),'999'))
        rejected(lambda:probe_native(job,descriptor,'/software-fixture',root/'proc'));negatives.append('changed_address_space_limit')
        (fake/'limits').write_text(limits)
    with patch('baliphy_sampler_resource_observation.psutil.Process',side_effect=psutil.NoSuchProcess(123)):
        assert probe_native(job,descriptor,'/software-fixture',root/'proc')['status']=='process_unavailable_during_observation'
    descriptors={};samples=[]
    for i,j in enumerate(jobs):
        cid=j['chain']['chain_id']
        # Include both missing live reads and missing native attempts. These
        # are artificial full-grid accounting cases, not measured production.
        if i==1619:continue
        descriptors[cid]=dict(chain_id=cid,pid=1000+i,created=100.,identity_sha256='software'+str(i))
        if i==1618:continue
        samples.append(dict(chain_id=cid,pid=1000+i,created=100.,identity_sha256='software'+str(i),
            checked_utc='2026-10-03T00:00:00+00:00',status='verified_live_native_resource_observation',
            reported_memory_bytes=expected,cpu_seconds=3.0))
    rows=finish_roles(jobs,descriptors,iter(samples),{})
    snapshots=[dict(cgroup=dict(memory_bytes={'memory.current':10,'memory.peak':20}),elapsed_since_previous_snapshot_seconds=1.)]
    summary=summarize(rows,iter(snapshots),'software_failure')
    assert summary['roles_with_native_attempt']==1619 and summary['roles_with_live_observation']==1618
    assert summary['roles_without_live_observation']==2 and summary['roles_without_native_attempt']==1
    assert summary['maximum_observed_cgroup_memory_bytes']==10 and summary['maximum_observed_cgroup_reported_peak_bytes']==20
    bad=dict(samples[0],pid=9999);rejected(lambda:finish_roles(jobs,descriptors,[bad],{}));negatives.append('reused_or_changed_sample_pid')
    bad=dict(samples[0],identity_sha256='changed');rejected(lambda:finish_roles(jobs,descriptors,[bad],{}));negatives.append('changed_sample_identity_hash')
    # One actual native software process checks Linux units, PID/creation,
    # argv after prlimit exec, exact caps, cgroup membership and live memory.
    prior=Path('metadata/baliphy_reference_sampler_software_validation_20261003.json')
    earlier=json.loads(prior.read_text());template=next(r for r in earlier['native_rows'] if r['prior_label']=='package')
    source_receipt=Path(template['native_receipt']);config=json.loads((source_receipt.parent.parent/'configuration.json').read_text())
    config['command'][config['command'].index('--seed')+1]='20263001'
    config['command'][config['command'].index('--cpu=300')]='--cpu=360'
    live_job=dict(chain=dict(chain_id='software-live-memory'),source_seed=425,
        memory_reservation_bytes=12*2**30,config=config)
    group=next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
    probes=[];folder=root/'native';identity=folder/'attempt-0001'/'process.json'
    with ThreadPoolExecutor(max_workers=1) as pool:
        future=pool.submit(run_attempt,folder,config)
        while not future.done():
            if identity.exists():probes.append(probe_native(live_job,describe(live_job,identity),group))
            time.sleep(.05)
        receipt=future.result()
    assert json.loads(receipt.read_text())['exit_code']==0
    live=[x for x in probes if x['status']=='verified_live_native_resource_observation'];assert live
    write_json(root/'live_native_probes.json',probes);write_json(root/'full_grid_accounting.json',summary)
    paths=[Path(__file__),Path('scripts/baliphy_sampler_resource_observation.py'),
        Path('scripts/run_baliphy_sampler_resource_observation.py'),native_path,Path(native['jobs']),prior,
        root/'live_native_probes.json',root/'full_grid_accounting.json',receipt,folder/'configuration.json']
    bindings={str(p):sha(p) for p in paths}
    for name,h in json.loads(receipt.read_text())['artifacts'].items():bindings[str(receipt.parent/name)]=h
    bindings.update(config['pins']);verify(dict(pins=bindings))
    result=dict(status='passed_full_sampler_resource_observation_software_contracts',
        checked_utc=datetime.now(timezone.utc).isoformat(),full_roles_checked=1620,
        full_missing_observation_accounting_checked=True,reused_pid_not_adopted=True,ephemeral_process_not_terminal_claim=True,
        wrapper_not_mistaken_for_native=True,actual_live_native_observations=len(live),
        empty_or_transitioning_argv_not_accepted_as_native=True,stable_unrelated_command_rejected=True,
        actual_native_receipt=str(receipt),malformed_or_altered_observations_rejected=negatives,
        source_hashes=bindings,scientific_eligibility=False,
        scope='Full1620real role metadata with artificial missing observations; one actual capped synthetic five-tip20iteration native process, sampled only while identity/caps/cgroup match. No fungal pilot, native production restart, exact final peak, long-chain memory qualification or scientific acceptance.')
    with args.receipt.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)


if __name__=='__main__':main()
