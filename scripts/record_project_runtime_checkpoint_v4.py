#!/usr/bin/env python3
"""Refresh original identities, preserved failures and full closed artifact/journal proofs."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess

import psutil

from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def fingerprint(record):
    try:
        proc=psutil.Process(record['pid'])
        if proc.create_time()!=record['created'] or proc.status()==psutil.STATUS_ZOMBIE:return None
        assert proc.cmdline()==record['cmdline'],record['pid']
        return proc
    except psutil.NoSuchProcess:return None


def journal_terminal(record,expected_failure=False):
    assert fingerprint(record) is None
    launch=json.loads(Path(record['launch']).read_text())
    unit=launch.get('unit',launch.get('service'))
    raw=subprocess.check_output(['journalctl','--user','-u',unit,'-o','json','--no-pager'],text=True)
    rows=[json.loads(line) for line in raw.splitlines()]
    exact=[r for r in rows if r.get('_PID')==str(record['pid']) and
           r.get('_CMDLINE')==' '.join(record['cmdline'])]
    assert exact,unit
    invocations={r['_SYSTEMD_INVOCATION_ID'] for r in exact}
    resources=[r for r in rows if r.get('USER_INVOCATION_ID') in invocations and r.get('CPU_USAGE_NSEC')]
    assert resources or expected_failure,unit
    failed=any(('Main process exited' in (r.get('MESSAGE') or '') or
                'Failed with result' in (r.get('MESSAGE') or ''))
                for r in rows if r.get('USER_INVOCATION_ID') in invocations)
    assert failed is expected_failure,unit
    return dict(**{k:record[k] for k in ['launch','pid','created','cmdline']},
        launch_sha256=sha(record['launch']),
        status=('verified_original_failure_preserved' if failed else
                'verified_original_terminal_success_with_bound_completed_artifacts'),
        journal_sha256_at_observation=hashlib.sha256(raw.encode()).hexdigest(),
        captured_process_messages=len(exact),
        original_failure_records=[{k:r.get(k) for k in ['__REALTIME_TIMESTAMP','MESSAGE','EXIT_CODE','EXIT_STATUS','RESULT','USER_INVOCATION_ID']}
            for r in rows if expected_failure and r.get('USER_INVOCATION_ID') in invocations and
            ('Main process exited' in (r.get('MESSAGE') or '') or 'Failed with result' in (r.get('MESSAGE') or ''))],
        completion_resource_records=[{
            k:r.get(k) for k in ['__REALTIME_TIMESTAMP','CPU_USAGE_NSEC','MEMORY_PEAK','MEMORY_SWAP_PEAK']}
            for r in resources])


def live_record(record,proc):
    result=dict(record)
    children=[]
    vanished=[]
    for child in proc.children(recursive=True):
        try:
            children.append(dict(pid=child.pid,created=child.create_time(),
                cmdline=child.cmdline(),status=child.status(),
                cpu_seconds_at_observation=sum(child.cpu_times()[:2])))
        except (psutil.NoSuchProcess,psutil.ZombieProcess) as error:
            vanished.append(dict(pid=child.pid,error_type=type(error).__name__,
                scope='Ephemeral child disappeared during observation; not a parent-job failure.'))
    result.update(status=proc.status(),cpu_seconds_at_observation=sum(proc.cpu_times()[:2]),
        children=children,children_unavailable_during_observation=vanished)
    if 'launch' in result:result['launch_sha256']=sha(result['launch'])
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--add-launch',action='append',default=[])
    parser.add_argument('--completion',action='append',default=[])
    parser.add_argument('--verified-receipt',action='append',default=[])
    parser.add_argument('--failed-launch',action='append',default=[])
    parser.add_argument('--coalescent-completion',default='metadata/species_coalescent_quartets_completed_20261002_v2.json')
    args=parser.parse_args()
    previous=json.loads(args.previous.read_text());bindings={}
    closed=dict(previous['full_closed_stages_verified'])
    completion_paths=[r['path'] for r in closed.values()]+args.completion
    for path in completion_paths:
        r=json.loads(Path(path).read_text())
        archive_path=r['full_hash_archive']
        bind(bindings,path);bind(bindings,archive_path,r['full_hash_archive_sha256'])
        proof=json.loads(Path(archive_path).read_text())
        for p,d in proof['source_hashes'].items():bind(bindings,p,d)
        assert len(proof['services'])==r['exact_process_journals_checked']
        closed[Path(path).stem]=dict(path=path,sha256=sha(path),status=r['status'],
            full_hash_archive=archive_path,full_hash_archive_sha256=r['full_hash_archive_sha256'],
            full_bindings_checked=len(proof['source_hashes']),
            original_completion_journals=len(proof['services']),
            summary={k:r[k] for k in ['cases','marker_states','original_gene_split_support_decisions',
                'internal_branches','branch_gene_states','global_numerator_independently_verified'] if k in r})
    for path in args.verified_receipt:
        r=json.loads(Path(path).read_text());assert r['status'].startswith('passed_')
        bind(bindings,path)
        for p,d in r['source_hashes'].items():bind(bindings,p,d)
    verify(bindings)
    current=dict(previous)
    current.update(status='verified_live_and_completed_project_runtime_checkpoint',
        checked_utc=datetime.now(timezone.utc).isoformat(),previous_checkpoint=str(args.previous),
        previous_checkpoint_sha256=sha(args.previous),full_closed_stages_verified=closed,
        total_distinct_closed_bindings_checked=len(bindings))
    handles=previous['verified_pipeline_handles']+previous['verified_terminal_pipeline_handles']
    for path in args.add_launch:
        r=json.loads(Path(path).read_text());r['launch']=path;handles.append(r)
    assert len({r['launch'] for r in handles})==len(handles)
    live,terminal,new_failures=[],[],[]
    failed_paths=set(args.failed_launch)
    assert failed_paths<={r['launch'] for r in handles}
    for r in handles:
        if 'launch' in r:
            launch=json.loads(Path(r['launch']).read_text())
            assert all(r[k]==launch[k] for k in ['pid','created','cmdline'])
            assert sha(launch['plan'])==launch['plan_sha256']
        proc=fingerprint(r)
        if r['launch'] in failed_paths:
            assert proc is None,('Declared failure remains live',r['launch'])
            new_failures.append(journal_terminal(r,expected_failure=True))
            continue
        if proc is not None:live.append(live_record(r,proc))
        else:terminal.append(journal_terminal(r))
    originals=[]
    for r in previous['preserved_original_live_jobs']:
        proc=fingerprint(r)
        assert proc is not None,('Original scientific job needs terminal audit',r['pid'])
        originals.append(live_record(r,proc))
    failures=[journal_terminal(r,expected_failure=True) for r in previous['preserved_failed_original_handles']]+new_failures
    current.update(verified_pipeline_handles=live,verified_terminal_pipeline_handles=terminal,
        preserved_original_live_jobs=originals,preserved_failed_original_handles=failures)
    state_path=Path(previous['native_progress']['path']);state=json.loads(state_path.read_text())
    assert state['total']==298848 and 0<=state['completed']<=state['total']
    current['native_progress']=dict(path=str(state_path),sha256_at_observation=sha(state_path),
        **{k:state[k] for k in ['stage','completed','total','counts']})
    native_root=Path('results/phylogeny/full-species-coalescent-sensitivities-20261002-v1')
    progress=sorted(native_root.glob('progress_*.json'))
    current['coalescent_native_progress']=dict(completed_cases=len(progress),total_cases=30,
        latest_case=json.loads(progress[-1].read_text())['latest_case'] if progress else None,
        full_batch_complete=(native_root/'receipt.json').is_file(),
        full_independent_quartet_closure_complete=Path(args.coalescent_completion).exists())
    current.update(coalescent_inputs_complete=True,first_full_coalescent_candidate_quartets_verified=True,
        all_eight_scientific_aims_incomplete=True,gpu_inference_paused=True,
        scope='Exact original PID/create/CMD checked for every ongoing pipeline and six original scientific/retrieval jobs; completed handles require actual invocation-linked original completion/resource journals and full closed source/artifact hashes. Initial failed readers remain preserved. Full coalescent inputs and first complete native-candidate independent global/local quartet/numeric proof are closed; exact native30-case output availability, full numerical closure and corrected downstream progress are recorded separately. Preserved original numerical/dependency failures are checked by exact invocation journals, not silently marked successful. Native constraints/model/gene uncertainty, roots, reconciliation and all eight biological aims remain open. GPU prediction remains paused; other projects may use the shared GPUs. No restarts or mutation of existing jobs.')
    with args.output.open('x') as f:f.write(json.dumps(current,indent=2)+'\n')
    print(json.dumps(dict(status=current['status'],live_pipeline_handles=len(live),
        original_live_jobs=len(originals),verified_terminal_successes=len(terminal),
        preserved_failures=len(failures),distinct_closed_bindings_checked=len(bindings),
        coalescent_progress=current['coalescent_native_progress'],background_progress=current['native_progress']['completed'])),flush=True)


if __name__=='__main__':main()
