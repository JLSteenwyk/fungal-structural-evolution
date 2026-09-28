"""Restartable complete-grid descriptive residual diagnostics, retaining failures."""
import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import numpy as np
import pandas as pd
import psutil
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha, write_json
from evaluate_selected_matched_residuals import evaluate

FACTORS = {}


def initialize(spec):
    threadpool_limits(1)
    for tree, source in spec.items():
        assert sha(source['path']) == source['sha256']
        with np.load(source['path'], allow_pickle=False) as h:
            FACTORS[tree] = h['factor']


def saved_result(path, binding):
    digest = path.with_suffix('.sha256')
    if not path.exists() and not digest.exists():
        return None
    assert path.exists() and digest.exists(), 'Incomplete checkpoint pair requires review'
    assert sha(path) == digest.read_text().strip(), 'Checkpoint checksum mismatch'
    result = json.loads(path.read_text())
    assert result['binding'] == binding, 'Checkpoint source/plan mismatch'
    assert set(result['fits']) == set(binding['selected_sources'])
    for tree, fit in result['fits'].items():
        assert fit['tree'] == tree and fit['fit_input_id'] == binding['fit_input_id']
        assert fit['selected_source_sha256'] == binding['selected_sources'][tree]['sha256']
    return result


def execute(task, rows, cache, output, plan_hash):
    identifier = task['fit_input_id']
    binding = dict(plan_sha256=plan_hash, **task)
    target = Path(output)/'inputs'/identifier[:2]/(identifier+'.json')
    assert set(rows) == set(task['selected_sources']) == set(FACTORS)
    cache_path = Path(cache)/task['cache_entry']['path']
    assert cache_path.resolve().is_relative_to(Path(cache).resolve())
    assert sha(cache_path) == task['cache_entry']['sha256']
    for tree, source in task['selected_sources'].items():
        assert sha(source['path']) == source['sha256'] == rows[tree]['selected_source_sha256']
    previous = saved_result(target, binding)
    if previous is not None:
        return {'fit_input_id': identifier, 'path': str(target), 'sha256': sha(target)}
    entry = task['cache_entry']
    path = Path(cache)/entry['path']
    assert path.resolve().is_relative_to(Path(cache).resolve())
    assert sha(path) == entry['sha256']
    with np.load(path, allow_pickle=False) as h:
        arrays = {k: h[k] for k in h.files}
    for key, field in [('matrix','values_sha256'), ('row_identity','ordered_identity_sha256')]:
        assert hashlib.sha256(arrays[key].tobytes()).hexdigest() == entry['recipe'][field]
    assert set(rows) == set(task['selected_sources']) == set(FACTORS)
    fits = {}
    start = time.monotonic()
    for tree, row in rows.items():
        source = task['selected_sources'][tree]
        # Integrity failures stop the job; numerical failures remain explicit dispositions.
        assert sha(source['path']) == source['sha256'] == row['selected_source_sha256']
        assert row['fit_input_id'] == identifier and row['tree'] == tree
        assert int(row['records']) == entry['records']
        try:
            fits[tree] = evaluate(row, arrays, FACTORS[tree])
        except (AssertionError, ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
            fits[tree] = dict(fit_input_id=identifier, tree=tree,
                selected_source_sha256=source['sha256'], selection=row['selection'],
                status='residual_computation_requires_review',
                error_type=type(error).__name__, error=str(error))
    result = dict(binding=binding, fits=fits, elapsed_seconds=time.monotonic()-start)
    target.parent.mkdir(parents=True, exist_ok=True)
    write_json(target, result)
    target.with_suffix('.sha256').write_text(sha(target)+'\n')
    saved_result(target, binding)
    return {'fit_input_id': identifier, 'path': str(target), 'sha256': sha(target)}


def wait_for_replay(spec):
    launch = json.loads(Path(spec['launch']).read_text())
    while True:
        state = dict(line.split('=',1) for line in subprocess.check_output(
            ['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','MainPID',
             '-p','Result','-p','ExecMainStatus'], text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:
            assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0'
            return
        assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid']
        try:
            process = psutil.Process(launch['pid'])
            assert process.create_time()==launch['created'] and process.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:
            time.sleep(1)
            continue
        print('Waiting for complete cache replay', flush=True)
        time.sleep(30)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--plan', type=Path, required=True)
    args=ap.parse_args(); config=json.loads(args.plan.read_text()); ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for name,digest in config['pins'].items():
            assert sha(name)==digest, name
    verify()
    root=Path(config['inputs']); preparation=json.loads((root/'receipt.json').read_text())
    wait_for_replay(preparation['required_prerequisite']); verify()
    replay_path=Path(preparation['required_prerequisite']['receipt'])
    replay=json.loads(replay_path.read_text())
    assert replay['status']==preparation['required_prerequisite']['status']
    assert replay['fits']==144040 and replay['inputs']==28808
    assert replay['original_receipt_sha256']==sha(config['original_receipt'])
    assert replay['plan_sha256']==sha(config['replay_plan'])
    assert replay['source_cache_receipt_sha256']==sha(Path(config['cache'])/'receipt.json')
    for name,digest in replay['artifacts'].items(): assert sha(replay_path.parent/name)==digest
    for name,digest in preparation['source_bindings'].items(): assert sha(name)==digest
    for name,digest in preparation['artifacts'].items(): assert sha(root/name)==digest
    tasks=[json.loads(line) for line in (root/'tasks.jsonl').open()]
    assert len(tasks)==28808 and len({t['fit_input_id'] for t in tasks})==28808
    selected=pd.read_parquet(config['selected'])
    assert len(selected)==144040 and not selected.duplicated(['fit_input_id','tree']).any()
    rows={}
    for row in selected.to_dict('records'):
        rows.setdefault(row['fit_input_id'], {})[row['tree']]=row
    assert set(rows)=={t['fit_input_id'] for t in tasks}
    out=Path(config['output']); out.mkdir(parents=True, exist_ok=True)
    lock=(out/'run.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
    stored_plan=out/'run_plan.json'
    if stored_plan.exists(): assert sha(stored_plan)==ph
    else: shutil.copyfile(args.plan, stored_plan)
    manifest=[]; iterator=iter(tasks); resources=config['resources']
    def disk_check():
        assert shutil.disk_usage(out).free>resources['free_disk_reserve_gib']*2**30
    disk_check()
    with ProcessPoolExecutor(max_workers=resources['workers'], initializer=initialize,
                             initargs=(preparation['factors'],)) as pool:
        pending=set()
        def submit():
            task=next(iterator,None)
            if task is None: return False
            pending.add(pool.submit(execute,task,rows[task['fit_input_id']],config['cache'],str(out),ph))
            return True
        for _ in range(resources['maximum_inflight_tasks']):
            if not submit(): break
        while pending:
            finished,_=wait(pending, return_when=FIRST_COMPLETED)
            for future in finished:
                pending.remove(future); manifest.append(future.result()); submit()
                if len(manifest)%100==0:
                    disk_check(); print('Completed descriptive residual inputs',len(manifest),'/28808',flush=True)
    assert len(manifest)==28808 and {r['fit_input_id'] for r in manifest}==set(rows)
    verify(); manifest.sort(key=lambda r:r['fit_input_id'])
    (out/'manifest.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in manifest))
    write_json(out/'receipt.json',dict(status='complete_descriptive_residual_dispositions_pending_full_readback',
        inputs=28808,fits=144040,plan_sha256=ph,replay_receipt_sha256=sha(replay_path),
        artifacts={'manifest.jsonl':sha(out/'manifest.jsonl'),'run_plan.json':sha(stored_plan)},
        scope='All selected fits attempted; descriptive marginal residual diagnostics with failures retained. '
              'Full output readback and simulation-based reference distributions remain required; no adequacy test claimed.'))


if __name__=='__main__': main()
