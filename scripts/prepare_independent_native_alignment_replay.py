#!/usr/bin/env python3
"""Replay every selected native chain with immutable whole-chain checkpoints."""
import argparse
import fcntl
import json
from pathlib import Path
import shutil
import time

from independent_native_alignment_replay_sources import load, summarize
from independent_native_alignment_replay import replay_chain
from prepare_independent_baliphy_scalar_readback_v2 import atomic
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def marker(path, source):
    return dict(schema='independent-native-alignment-replay-v1', plan_sha256=sha(path),
        inventory_details_sha256=source['inventory']['full_chain_details_sha256'])


def artifacts(root, entries, results):
    paths = [root/'stage_plan.json',root/'chain_manifest.json',*(root/e['path'] for e in entries)]
    paths.extend(root/'chains'/row['frames_path'] for row in results.values() if 'frames_path' in row)
    assert len(paths) == len(set(paths))
    assert set((root/'chains').glob('*.json')) == {root/e['path'] for e in entries}
    assert set((root/'chains').glob('*.jsonl')) == {p for p in paths if p.suffix=='.jsonl'}
    return {str(p.relative_to(root)):sha(p) for p in paths}


def run(path, stop_after_chains=None):
    started=time.perf_counter();plan=json.loads(path.read_text());root=Path(plan['output'])
    assert not (root/'receipt.json').exists(), 'Completed producer cannot restart'
    assert shutil.disk_usage(root.parent).free >= plan['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(exist_ok=True); lock=(root/'stage.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (root/'receipt.json').exists(), 'Completed producer cannot restart'
    source,bindings=load(plan,path);stage=marker(path,source);mp=root/'stage_plan.json'
    if mp.exists(): assert json.loads(mp.read_text())==stage
    else: atomic(mp,stage)
    folder=root/'chains';folder.mkdir(exist_ok=True);results={};entries=[]
    for number,cid in enumerate(sorted(source['details']),1):
        target=folder/(cid+'.jsonl');cp=folder/(cid+'.json')
        expected=replay_chain(source,cid,target,readback=target.exists())
        document=dict(stage=stage,chain_id=cid,result=expected)
        if cp.exists():assert json.loads(cp.read_text())==document
        else:atomic(cp,document)
        results[cid]=expected;entries.append(dict(chain_id=cid,path=str(cp.relative_to(root)),sha256=sha(cp)))
        print('independent_native_alignment_chains',number,'/',len(source['details']),
            'seconds',round(time.perf_counter()-started,3),flush=True)
        if number==stop_after_chains:raise InterruptedError('Software checkpoint contract')
    summary=summarize(results,source,plan);atomic(root/'chain_manifest.json',entries)
    verify(bindings);exported=artifacts(root,entries,results)
    result=dict(status='complete_full_independent_native_alignment_replay_pending_readback',
        plan_sha256=sha(path),**summary,source_hashes=bindings,artifacts=exported,
        elapsed_seconds=time.perf_counter()-started,scientific_eligibility=False,scope=plan['scope'])
    atomic(root/'receipt.json',result);print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();run(a.plan)
