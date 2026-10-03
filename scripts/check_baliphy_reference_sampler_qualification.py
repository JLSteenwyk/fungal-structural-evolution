#!/usr/bin/env python3
"""Full resource scheduler contracts and native synthetic five-tip output checks."""
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import threading
import time
from unittest.mock import patch

from ancestral_chain_attempt import sha, write_json
from baliphy_reference_initialization import transform
from baliphy_reference_sampler_qualification import inspect, summarize
from prepare_baliphy_reference_sampler_qualification import inputs
from reference_sampler_memory_budget import MemoryBudget
import run_baliphy_reference_sampler_qualification as workflow


def rejected(action):
    try:action()
    except (AssertionError,KeyError,ValueError,RuntimeError):return
    raise AssertionError('Altered or unsafe qualification accepted')


def wait_for_waiter(budget, identifier):
    end=time.monotonic()+5
    while time.monotonic()<end:
        with budget.condition:
            if identifier in budget.waiting:return
        time.sleep(.005)
    raise AssertionError('Expected reservation waiter did not arrive')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--receipt',type=Path,required=True);args=parser.parse_args()
    root=args.output.resolve();root.mkdir(exist_ok=False)
    startup,jobs,risk,source_paths=inputs(); negatives=[]
    assert len(jobs)==1620 and risk==['OG0000972','OG0001082']
    assert sum(j['memory_reservation_bytes']==48*2**30 for j in jobs)==108
    for job in jobs:
        cmd=job['config']['command']
        assert '--test' not in cmd and cmd[cmd.index('--iterations')+1]=='20'
        assert '--as='+str(job['memory_reservation_bytes']) in cmd
        assert '--cpu='+str(job['config']['timeout_seconds']) in cmd
        assert job['config']['pins'][job['chain']['program']]==sha(job['chain']['program'])
    budget=MemoryBudget(192*2**30)
    def software_role(job):
        with budget.reserve(job['chain']['chain_id'],job['memory_reservation_bytes']):time.sleep(.001)
    with ThreadPoolExecutor(max_workers=16) as pool:list(pool.map(software_role,jobs))
    stress=workflow.verify_ledger(budget.events,jobs,192*2**30,16)
    write_json(root/'full_scheduler_stress_events.json',budget.events)
    for name in ['missing_release','duplicate_acquire','changed_reservation','oversubscribed_total']:
        events=copy.deepcopy(budget.events)
        if name=='missing_release':events.pop()
        elif name=='duplicate_acquire':events[1]=copy.deepcopy(events[0])
        elif name=='changed_reservation':events[0]['amount']+=1
        else:events[0]['reserved_total']=193*2**30
        rejected(lambda:workflow.verify_ledger(events,jobs,192*2**30,16));negatives.append(name)
    fifo=MemoryBudget(192);holders=[fifo.reserve('held'+str(i),x) for i,x in enumerate([48,48,48,12])]
    for hold in holders:hold.__enter__()
    heavy_entered=threading.Event();heavy_release=threading.Event();light_entered=threading.Event()
    def heavy():
        with fifo.reserve('heavy',48):heavy_entered.set();assert heavy_release.wait(5)
    def light():
        with fifo.reserve('light',12):light_entered.set()
    with ThreadPoolExecutor(max_workers=2) as pool:
        hf=pool.submit(heavy);wait_for_waiter(fifo,'heavy');lf=pool.submit(light);wait_for_waiter(fifo,'light')
        assert not heavy_entered.is_set() and not light_entered.is_set()
        holders[3].__exit__(None,None,None);assert heavy_entered.wait(5);assert not light_entered.is_set()
        holders[1].__exit__(None,None,None);assert light_entered.wait(5);heavy_release.set();hf.result();lf.result()
    for i in [0,2]:holders[i].__exit__(None,None,None)
    assert fifo.used==0
    abort=MemoryBudget(48);lease=abort.reserve('active',48);lease.__enter__()
    with ThreadPoolExecutor(max_workers=1) as pool:
        def pending():
            with abort.reserve('waiting',48):raise AssertionError('Aborted waiter admitted')
        future=pool.submit(pending);wait_for_waiter(abort,'waiting');abort.abort('Synthetic orphan/controller failure')
        rejected(future.result);lease.__exit__(None,None,None)
    assert abort.used==0 and not abort.active
    rejected(lambda:abort.reserve('later',12).__enter__());negatives.append('new_admission_after_abort')
    recovery=MemoryBudget(12)
    def raising():
        with recovery.reserve('raises',12):raise ValueError('Synthetic audit exception')
    rejected(raising);assert recovery.used==0
    with recovery.reserve('later',12):assert recovery.used==12
    alignment=root/'alignment.faa';tree=root/'tree.nwk';mapping=root/'mapping.tsv'
    alignment.write_text('>a\nA-CDE-FG\n>b\nA-C-EYFG\n>c\nATCDE-F-\n>d\n--CDEF-G\n>e\nAT-DE-FG\n')
    tree.write_text('((((a:0.1,b:0.1):0.1,c:0.1):0.1,d:0.1):0.1,e:0.1);\n')
    with mapping.open('w') as f:
        writer=csv.DictWriter(f,delimiter='\t',fieldnames=['guide','family','dataset','level','source_node','retained_set_json'])
        writer.writeheader()
        for i,n in enumerate([2,3,4,5]):writer.writerow(dict(guide='profile',family='software',dataset='whole',level=str(i),
            source_node='software'+str(i),retained_set_json=json.dumps(list('abcde')[:n])))
    native_rows=[];fixtures=root/'native';fixtures.mkdir();(fixtures/'chains').mkdir();bindings={}
    native_budget=MemoryBudget(12*2**30)
    for prior in ['broad','centered','package']:
        original=next(x['chain'] for x in json.loads(Path(startup['jobs']).read_text()) if x['chain']['prior_label']==prior)
        source=Path(original['program']).read_text().replace(str(Path(original['alignment']).resolve()),str(alignment)).replace(str(Path(original['tree']).resolve()),str(tree))
        program=root/(prior+'.hs');program.write_text(transform(source))
        binary=Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/bali-phy').resolve()
        config=dict(command=['/usr/bin/prlimit','--as='+str(12*2**30),'--cpu=300','--fsize='+str(2*2**30),
            '--',str(binary),'--seed','426','run',str(program),'--iterations','20','--log-format','json,tsv','--name','independent-chain'],
            timeout_seconds=360,pins={str(x):sha(x) for x in [Path('/usr/bin/prlimit'),binary,program,alignment,tree]})
        chain=dict(original,chain_id='software-'+prior,seed=426,family='software',proteins=5,
            original_configuration_ids=['software-whole-'+prior],alignment=str(alignment),alignment_sha256=sha(alignment),
            tree=str(tree),tree_sha256=sha(tree),program=str(program),program_sha256=sha(program))
        job=dict(chain=chain,config=config,source_seed=425,memory_reservation_bytes=12*2**30)
        # Separate per-prior budgets permit this software fixture seed, without
        # changing the production globally distinct role seeds.
        row=workflow.execute_job(job,fixtures,'software-only',mapping,native_budget,1)
        assert row['status']=='full_short_sampler_output_integrity_checked_not_posterior'
        native_rows.append(row);receipt=Path(row['native_receipt']);bindings[str(receipt)]=sha(receipt)
        bindings[str(receipt.parent.parent/'configuration.json')]=sha(receipt.parent.parent/'configuration.json')
        for name,h in json.loads(receipt.read_text())['artifacts'].items():bindings[str(receipt.parent/name)]=h
        assert inspect(job,receipt,'software-only',mapping)==row
        print('native_five_tip_sampler_fixture',prior,'saved',row['saved_alignments'],'candidates',row['candidate_frames'],flush=True)
        # Actual native files/receipts stay immutable. Failed partial traces are
        # explicit mocked decoder cases, not invented native crash evidence.
        fake=json.loads(receipt.read_text());fake['exit_code']=1;fake['status']='failed'
        actual_loads=json.loads
        def patched_loads(text,*a,**kw):return fake if text==receipt.read_text() else actual_loads(text,*a,**kw)
        for name,error in [('failed_header_only_trace',AssertionError('No logged rows')),
                           ('failed_malformed_trace',ArithmeticError('Malformed scalar token'))]:
            with patch('baliphy_reference_sampler_qualification.json.loads',side_effect=patched_loads),patch(
                'baliphy_reference_sampler_qualification.trace_summary',side_effect=error):
                failed=inspect(job,receipt,'software-only',mapping)
                assert failed['status']=='unsuccessful_sampler_qualification_attempt_retained'
                assert failed['partial_scalar_trace'] is None and failed['partial_scalar_trace_issue']
        # No repeated native attempt, including after an existing failed receipt.
        with patch.object(workflow,'run_attempt',side_effect=AssertionError('Unexpected automatic retry')):
            fresh_budget=MemoryBudget(12*2**30)
            assert workflow.execute_job(job,fixtures,'software-only',mapping,fresh_budget,1)==row
    fake_rows=[]
    for i,job in enumerate(jobs):
        c=job['chain'];failed=i in [0,4]
        fake_rows.append(dict(chain_id=c['chain_id'],effective_input_group=c['effective_input_group'],
            model_input_identity=c['effective_input_group']+'-'+c['prior_label'],chain_role=c['chain'],
            original_configuration_ids=c['original_configuration_ids'],
            status='unsuccessful_sampler_qualification_attempt_retained' if failed else 'full_short_sampler_output_integrity_checked_not_posterior',
            partial_scalar_trace=None,saved_alignments=0 if failed else 3,candidate_frames=0 if failed else 12,
            allocation_warning_lines=0,bad_alloc=failed,posterior_qualified=False))
    summary=summarize(fake_rows)
    assert (summary['checked_sampler_attempts'],summary['unsuccessful_sampler_attempts'],
        summary['complete_quartets'],summary['unresolved_quartets'])==(1618,2,403,2)
    for name,change in [('missing_role',lambda x:x.pop()),('duplicate_role',lambda x:x.__setitem__(1,x[0])),
                        ('changed_chain_role',lambda x:x[0].update(chain_role=2))]:
        bad=copy.deepcopy(fake_rows);change(bad);rejected(lambda:summarize(bad));negatives.append(name)
    sources=[Path(__file__),Path('scripts/reference_sampler_memory_budget.py'),
        Path('scripts/baliphy_reference_sampler_qualification.py'),Path('scripts/run_baliphy_reference_sampler_qualification.py'),
        Path('scripts/prepare_baliphy_reference_sampler_qualification.py'),Path('scripts/readback_independent_baliphy_chain.py'),
        Path('scripts/ancestral_chain_attempt.py'),alignment,tree,mapping,*source_paths]
    bindings.update({str(p):sha(p) for p in sources})
    result=dict(status='passed_full_reference_sampler_qualification_software_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        full_chain_configurations_checked=1620,stress_roles=1620,stress_reservation_audit=stress,
        fifo_large_role_not_bypassed=True,admission_abort_blocks_waiters=True,exception_releases_reservation=True,
        native_fixture_priors=['broad','centered','package'],native_saved_alignments_checked=9,native_candidate_frames_checked=36,
        failed_header_only_trace_retained=True,failed_malformed_trace_retained=True,
        simulated_failed_roles_retained=2,simulated_unresolved_quartets=2,
        altered_resource_ledgers_and_grid_accounting_rejected=negatives,native_rows=native_rows,
        source_hashes=bindings,scientific_eligibility=False,
        scope='Complete1620real production configurations/resource scheduler stress, FIFO/no-oversubscription/abort contracts. Three capped native synthetic five-tip20iteration fixtures test exact attempt controller and scalar/tree/tip/saved-alignment/four ancestral-candidate readback. Failed partial traces deliberately mocked; full metadata failure accounting synthetic. No fungal pilot, allocation repair, longer-chain or posterior qualification.')
    with args.receipt.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','native_rows']}),flush=True)


if __name__=='__main__':main()
