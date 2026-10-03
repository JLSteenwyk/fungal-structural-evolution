#!/usr/bin/env python3
"""Rebuild all measured resource/numerical/proposed-grid records from sources."""
import argparse
import fcntl
import json
from pathlib import Path

from ancestral_chain_attempt import write_json
from baliphy_horizon_resource_inventory import collect
from run_ortholog_pair_guide_comparison import sha


def verify_exports(root,data):
    root=Path(root)
    for name in ['chains','attempts','quartets','proposed']:
        with (root/(name+'.jsonl')).open() as f:rows=[json.loads(line) for line in f]
        assert rows==data[name],name


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());root=Path(plan['output']);digest=sha(a.plan)
    lock=(root/'readback.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (root/'readback.json').exists(),'Completed readback cannot be rerun'
    receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='complete_full_baliphy_horizon_resource_inventory_pending_readback' and receipt['plan_sha256']==digest
    for name,d in receipt['artifacts'].items():assert sha(root/name)==d
    data,bindings=collect(plan);bindings[str(a.plan)]=digest
    assert bindings==receipt['source_hashes']
    verify_exports(root,data)
    assert all(receipt[k]==v for k,v in data['summary'].items())
    for path,d in bindings.items():assert sha(path)==d,path
    for name,d in receipt['artifacts'].items():assert sha(root/name)==d
    write_json(root/'readback.json',dict(status='passed_full_baliphy_horizon_resource_inventory_serialized_readback',
        **data['summary'],plan_sha256=digest,producer_receipt=str(root/'receipt.json'),producer_receipt_sha256=sha(root/'receipt.json'),
        source_hashes=bindings,artifacts=receipt['artifacts'],
        scope='Full source reconstruction of all1620chains/1623attempts/405quartets/fresh-seed proposals and serialized records. '
        'Uses the same new census implementation, not a third algorithm. Native samples are sized only, not freshly decoded. '
        'No native sampling/GPU/restarts/qualified ensemble; all eight aims incomplete.'))
    print(json.dumps(data['summary'],allow_nan=False))


if __name__=='__main__':main()
