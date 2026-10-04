#!/usr/bin/env python3
"""Compare all24 failed roles once; preserve native failures and independently reread every outcome."""
import argparse
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import json
from pathlib import Path

import psutil

from ancestral_chain_attempt import sha,write_json
from baliphy_joint_fasta_v7_failure_jobs_v1 import failed_sources,validate_jobs
from baliphy_joint_sampler_scalar_v6 import inspect,SUCCESS
from reference_measurement_union_sources import bind,verify
from reference_sampler_memory_budget import MemoryBudget
from run_baliphy_scalar_v6_sampler_v1 import execute_job,verify_ledger,runtime_caps


def gates(plan):
    pins={}
    for gate in plan['gates']:
        receipt=Path(gate['receipt']);transport=Path(gate['transport'])
        r,t=[json.loads(p.read_text()) for p in [receipt,transport]]
        assert r['status']==gate['status']
        assert t['validation_sha256']==sha(receipt) and t['original_tool_terminal_exit_code']==0
        assert t['whole_wrapper_initial_and_terminal_payloads_matched']
        verify(t['source_hashes'])
        for path,h in t['source_hashes'].items():bind(pins,path,h)
        for path in [receipt,transport]:bind(pins,path)
    return pins


def summarize(rows):
    assert len(rows)==len({r['chain_id'] for r in rows})==24
    assert all(r['scientific_eligibility'] is r['posterior_qualified'] is False for r in rows)
    groups=defaultdict(list)
    for row in rows:groups[(row['effective_input_group'],row['prior_label'])].append(row)
    assert len(groups)==6 and all(sorted(r['chain_role'] for r in group)==[1,2,3,4] for group in groups.values())
    frames=[frame for row in rows for frame in row['joint_frames']]
    return dict(roles=24,prior_quartets=6,effective_inputs=2,
        status_counts=dict(Counter(r['status'] for r in rows)),
        native_zero_exit_roles=sum(r['exit_code']==0 for r in rows),
        integrity_checked_roles=sum(r['status']==SUCCESS for r in rows),
        complete_integrity_quartets=sum(all(r['status']==SUCCESS for r in group) for group in groups.values()),
        scalar_mapped_values=sum((r['scalar_v6_audit'] or {}).get('mapped_values_compared',0) for r in rows),
        saved_joint_frames=len(frames),ancestral_residue_category_pairs=sum(f['ancestral_pairs'] for f in frames),
        tip_residue_category_pairs=sum(f['tip_pairs'] for f in frames),
        full_original_failure_grid_output_integrity_checked=all(r['status']==SUCCESS for r in rows),
        posterior_qualified=False)


def run(path,reader=False):
    plan=json.loads(path.read_text());verify(plan['pins']);digest=sha(path)
    original,pins,mapping=failed_sources()
    jobs=json.loads(Path(plan['jobs']).read_text());scope=validate_jobs(jobs,original)
    assert scope==plan['job_scope'] and mapping==plan['mapping']
    resources=plan['resources']
    assert resources['cpus']==resources['workers']==4
    assert resources['reservation_capacity_gib']==192 and resources['memory_gib']==200
    limits=runtime_caps(plan,reader)
    pins.update(plan['pins']);pins.update(gates(plan))
    root=Path(plan['output']);started=datetime.now(timezone.utc).isoformat()
    process=psutil.Process()
    if reader:
        assert process.rlimit(psutil.RLIMIT_AS)==(24*2**30,)*2
        process.rlimit(psutil.RLIMIT_AS,(8*2**30,24*2**30))
        producer_path=root/'receipt.json';producer=json.loads(producer_path.read_text())
        assert producer['status']=='complete_all24_v7_failure_comparison_dispositions_pending_readback'
        assert producer['plan_sha256']==digest and not (root/'readback.json').exists()
        verify(producer['source_hashes'])
        for name,h in producer['artifacts'].items():assert sha(root/name)==h
        rows=[]
        for job in jobs:
            cid=job['chain']['chain_id'];cp=root/'chains'/(cid+'.json');saved=json.loads(cp.read_text())
            row=inspect(job,Path(saved['native_receipt']),digest,mapping,root/'frames'/cid,allow_export_creation=False)
            assert row==saved;rows.append(row)
            print('v7_failure_grid_readback',len(rows),'/24',cid,row['status'],flush=True)
        rows.sort(key=lambda r:r['chain_id'])
        assert rows==json.loads((root/'dispositions.json').read_text())
        assert plan==json.loads((root/'stage_plan.json').read_text())
    else:
        assert process.rlimit(psutil.RLIMIT_AS)==(192*2**30,)*2
        process.rlimit(psutil.RLIMIT_AS,(8*2**30,192*2**30))
        assert psutil.virtual_memory().available>=resources['memory_gib']*2**30
        root.mkdir(exist_ok=False);(root/'chains').mkdir()
        write_json(root/'stage_plan.json',plan)
        journal=root/'memory_reservations.jsonl';rows=[]
        with journal.open('x') as handle:
            def persist(event):handle.write(json.dumps(event,allow_nan=False)+'\n');handle.flush()
            budget=MemoryBudget(192*2**30,persist)
            with ThreadPoolExecutor(max_workers=4) as pool:
                for future in as_completed([pool.submit(execute_job,j,root,digest,mapping,budget,
                                            resources['minimum_free_disk_gib']) for j in jobs]):
                    row=future.result();rows.append(row)
                    print('v7_failure_grid_outcome',len(rows),'/24',row['chain_id'],row['status'],flush=True)
            assert budget.aborted is None
        rows.sort(key=lambda r:r['chain_id']);write_json(root/'dispositions.json',rows)
    events=[json.loads(line) for line in (root/'memory_reservations.jsonl').read_text().splitlines()]
    ledger=verify_ledger(events,jobs,192*2**30,4);summary=summarize(rows)
    for row in rows:
        receipt=Path(row['native_receipt']);bind(pins,receipt);bind(pins,receipt.parent.parent/'configuration.json')
        for name,h in json.loads(receipt.read_text())['artifacts'].items():bind(pins,receipt.parent/name,h)
        for frame in row['joint_frames']:bind(pins,frame['projection_array'],frame['projection_array_sha256'])
    expected={Path(frame['projection_array']) for row in rows for frame in row['joint_frames']}
    assert set((root/'frames').glob('*/*'))==expected
    for p in [path,Path(__file__),root/'stage_plan.json',root/'dispositions.json',root/'memory_reservations.jsonl',
              *sorted((root/'chains').glob('*.json'))]:bind(pins,p)
    verify(pins);assert sha(path)==digest
    result=dict(status='passed_independent_all24_v7_failure_comparison_readback' if reader else
        'complete_all24_v7_failure_comparison_dispositions_pending_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(),started_utc=started,plan_sha256=digest,
        **summary,reservation_audit=ledger,actual_cgroup_limits=limits,
        controller_address_space_limits=process.rlimit(psutil.RLIMIT_AS),source_hashes=pins,
        scientific_eligibility=False,gpu=False,original_jobs_restarted=False,scope=plan['scope'])
    if reader:
        assert summary==producer['summary'] and ledger==producer['reservation_audit']
        result['producer_receipt_sha256']=sha(producer_path)
    else:
        result['summary']=summary
        result['artifacts']={str(p.relative_to(root)):sha(p) for p in [root/'stage_plan.json',root/'dispositions.json',
                            root/'memory_reservations.jsonl',*sorted((root/'chains').glob('*.json'))]}
    target=root/('readback.json' if reader else 'receipt.json')
    with target.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(summary),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--reader',action='store_true');a=p.parse_args();run(a.plan,a.reader)
