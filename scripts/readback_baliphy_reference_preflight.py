#!/usr/bin/env python3
"""Re-read every native startup/receipt/source and all full-grid dispositions."""
import argparse
import fcntl
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from run_baliphy_reference_preflight import inspect, summarize, verify


def readback(path):
    plan = json.loads(path.read_text()); digest = sha(path); verify(plan)
    root = Path(plan['output']); lock = (root / 'readback.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (root / 'readback.json').exists(), 'Completed readback cannot restart'
    producer_path = root / 'receipt.json'; producer = json.loads(producer_path.read_text())
    assert producer['status'] == 'complete_full_reference_startup_dispositions_pending_readback'
    assert producer['plan_sha256'] == digest
    bindings = dict(plan['pins']); bindings[str(path)] = digest; bindings[str(producer_path)] = sha(producer_path)
    for name, h in producer['artifacts'].items():
        assert sha(root / name) == h; bindings[str(root / name)] = h
    exported = json.loads((root / 'dispositions.json').read_text()); rows = []
    saved = {r['chain_id']: r for r in exported}; assert len(saved) == len(exported) == 1620
    jobs = json.loads(Path(plan['jobs']).read_text())
    assert saved.keys() == {j['chain']['chain_id'] for j in jobs}
    for number, job in enumerate(sorted(jobs, key=lambda x: x['chain']['chain_id']), 1):
        cid = job['chain']['chain_id']; cp = root / 'chains' / (cid + '.json')
        row = json.loads(cp.read_text()); assert row == saved[cid]
        receipt_path = Path(row['attempt_receipt']); receipt = json.loads(receipt_path.read_text())
        config_path = receipt_path.parent.parent / 'configuration.json'
        assert json.loads(config_path.read_text()) == job['config']
        actual = inspect(job, receipt_path, digest); assert actual == row, 'Serialized native startup differs'
        bindings[str(receipt_path)] = sha(receipt_path); bindings[str(config_path)] = sha(config_path)
        for name, h in receipt['artifacts'].items(): bindings[str(receipt_path.parent / name)] = h
        rows.append(actual)
        if number % 50 == 0: print('reference_startup_serialized_readback', number, '/1620', flush=True)
    summary = summarize(rows)
    assert all(producer[k] == v for k, v in summary.items())
    for name, h in bindings.items(): assert sha(name) == h, name
    result = dict(status='passed_full_reference_startup_serialized_readback', plan_sha256=digest,
        producer_receipt_sha256=sha(producer_path), **summary, source_hashes=bindings,
        scientific_eligibility=False, scope='Complete native stdout, runtime trees, sources, attempts, configurations, serialized checkpoints and full design reread. Native-vs-Python homology comparison and original degree-aware density checked again. Reader shares the startup decoder; not a third independent likelihood algorithm. ' + plan['scope'])
    with (root / 'readback.json').open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(summary), flush=True); return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args(); readback(a.plan)
