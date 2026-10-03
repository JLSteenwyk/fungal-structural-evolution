#!/usr/bin/env python3
"""Publish full measured initial-horizon census and an unlaunched longer grid."""
import argparse
import fcntl
import json
from pathlib import Path
import time

from ancestral_chain_attempt import write_json
from baliphy_horizon_resource_inventory import collect
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());digest=sha(a.plan)
    root=Path(plan['output']);root.mkdir(parents=True,exist_ok=False)
    lock=(root/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    started=time.monotonic();data,bindings=collect(plan);bindings[str(a.plan)]=digest
    for name in ['chains','attempts','quartets','proposed']:
        with (root/(name+'.jsonl')).open('x') as f:
            for row in data[name]:f.write(json.dumps(row,sort_keys=True,allow_nan=False)+'\n')
    assert sha(a.plan)==digest
    for path,d in bindings.items():assert sha(path)==d,path
    assert data['summary']['full_chains']==1620 and data['summary']['production_launch_allowed'] is False
    write_json(root/'receipt.json',dict(status='complete_full_baliphy_horizon_resource_inventory_pending_readback',
        plan_sha256=digest,**data['summary'],source_hashes=bindings,
        artifacts={p.name:sha(p) for p in root.glob('*.jsonl')},elapsed_seconds=time.monotonic()-started,
        scope='All1620original IDs/405quartets/135inputs and1623initial+selected-recovery attempts. '
        'Fresh small-source/scalar-log/stderr/configuration/receipt hashes and complete numeric-row survey. '
        'Large native alignment file sizes are observed without rehashing their contents; memory peaks unavailable. '
        'Longer horizon has all1620fresh seeds, retains both failures and aliases, and remains unlaunched. '
        'Linear runtime/storage projections are conditional descriptive scenarios, not forecasts, resource guarantees or convergence.'))
    print(json.dumps(data['summary'],allow_nan=False))


if __name__=='__main__':main()
