#!/usr/bin/env python3
"""Verify new uniform audit queue plus exact original upstream/recovery handles."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
from record_project_runtime_checkpoint_v4 import fingerprint,live_record,journal_terminal
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    previous_path=Path('metadata/full_entity_operator_queue_checkpoint_20261002_v4.json')
    previous=json.loads(previous_path.read_text());pins={}
    plan_path=Path('metadata/full_uniform_covariance_qualification_plan_20261002.json');plan=json.loads(plan_path.read_text())
    for path,d in plan['pins'].items():assert sha(path)==d; pins[path]=d
    prior=json.loads(Path('metadata/full_entity_operator_queue_checkpoint_20261002.json').read_text())
    records={r['launch']:r for r in prior['handles']}
    launches=json.loads(Path('metadata/full_uniform_covariance_qualification_launches_20261002.json').read_text())
    assert launches['source_plan_sha256']==sha(plan_path)
    for path in launches['launches']:
        r=json.loads(Path(path).read_text());r['launch']=path;assert path not in records;records[path]=r
    handles=[]
    for path,r in records.items():
        launch=json.loads(Path(path).read_text());assert all(r[k]==launch[k] for k in ['pid','created','cmdline'])
        assert sha(launch['plan'])==launch['plan_sha256']
        proc=fingerprint(r)
        if proc is not None:
            observed=live_record(r,proc)
            if path in launches['launches']:
                cg=next(line.split('::',1)[1] for line in Path('/proc',str(proc.pid),'cgroup').read_text().splitlines() if line.startswith('0::'))
                actual={name:Path('/sys/fs/cgroup',cg.lstrip('/'),name).read_text().strip()
                        for name in ['cpu.max','memory.max','memory.swap.max']}
                assert actual=={'cpu.max':'200000 100000','memory.max':str(32*2**30),'memory.swap.max':'0'}
                observed['actual_cgroup_limits']=actual
        else:
            observed=journal_terminal(r);observed['status']='verified_original_terminal_success_resource_journal'
        handles.append(observed)
    closed={}
    for name in ['model_inputs','operator_bank']:
        stage=previous['stages'][name];assert stage['accepted']
        assert sha(stage['path'])==stage['sha256']
        assert sha(stage['full_hash_archive'])==stage['full_hash_archive_sha256']
        closed[name]=stage
    result=dict(status='verified_full_uniform_covariance_original_queue_and_upstream_closures',
        checked_utc=datetime.now(timezone.utc).isoformat(),source_plan=str(plan_path),source_plan_sha256=sha(plan_path),
        prior_full_closure_checkpoint=str(previous_path),prior_full_closure_checkpoint_sha256=sha(previous_path),
        pins_verified=pins,handles=handles,sealed_completed_stages=closed,
        full_design_production_accepted=False,full_uniform_covariance_production_accepted=False,
        all_eight_scientific_aims_incomplete=True,gpu_inference_paused=True,
        source_hashes={str(Path(__file__)):sha(Path(__file__)),
            'scripts/record_project_runtime_checkpoint_v4.py':sha('scripts/record_project_runtime_checkpoint_v4.py')},
        scope='Exact original handles, immutable plans and source pins checked. New queue cgroup limits '
              'revalidated directly. Input/operator sealed closure archives retain their full-source '
              'hash proofs and actual original journals verified by the linked prior checkpoint. '
              'Design and uniform qualification acceptance remain pending. No restart or GPU use.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(status=result['status'],original_handles=len(handles),pins=len(pins),closed_stages=list(closed))))


if __name__=='__main__':main()
