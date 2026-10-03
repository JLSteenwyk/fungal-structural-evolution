#!/usr/bin/env python3
"""Run all 931 qualified matched-source fixed-topology point fits once."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
import fcntl
import json
from pathlib import Path
import shutil

from matched_predictor_branch_fits import load,execute,ROLES
from matched_predictor_branch_inputs import verify
from ancestral_chain_attempt import sha

STATUS='complete_full_matched_predictor_native_point_fits_pending_independent_readback'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());source,bindings=load(plan,a.plan);root=Path(plan['output'])
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(exist_ok=False);lock=(root/'stage.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    (root/'roles').mkdir();(root/'native').mkdir()
    jobs=[(key,cfg,role) for key,cfg in sorted(source['configs'].items()) for role in ROLES]
    assert len(jobs)==931;rows=[]
    with ThreadPoolExecutor(max_workers=plan['resources']['workers']) as pool:
        futures={pool.submit(execute,plan,source['root'],root,key,cfg,role,source['axes']['positions']):(key,role[0])
                 for key,cfg,role in jobs}
        for f in as_completed(futures):
            row=f.result();rows.append(row)
            assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
            print('matched_predictor_native_role',len(rows),'/',len(jobs),row['input_id'],row['role'],row['status'],flush=True)
    counts=Counter(r['status'] for r in rows);by_input={k:[] for k in source['configs']}
    for row in rows:by_input[row['input_id']].append(row)
    intact=sum(len(v)==7 and all(r['point_estimate'] is not None for r in v) for v in by_input.values())
    summary=dict(input_comparison_cases=source['completion']['comparison_cases'],unique_ready_inputs=len(source['configs']),
        native_roles=len(rows),native_status_counts=dict(counts),intact_inputs=intact,unresolved_inputs=len(source['configs'])-intact,
        point_branch_values=sum(len(r['point_estimate']['branches']) for r in rows if r['point_estimate'] is not None),
        native_seconds_sum=sum(r['elapsed_seconds'] for r in sorted(rows,key=lambda r:(r['input_id'],r['role']))))
    verify(bindings)
    artifacts={str(p.relative_to(root)):sha(p) for p in root.rglob('*') if p.is_file() and p.name!='stage.lock'}
    receipt=dict(status=STATUS,plan_sha256=sha(a.plan),**summary,source_hashes=bindings,artifacts=artifacts,
                 scientific_eligibility=False,scope=plan['scope'])
    with (root/'receipt.json').open('x') as f:f.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
