#!/usr/bin/env python3
"""Observe original V6 startup controllers, native limits and unclosed progress."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists()
    plan_path = Path('metadata/baliphy_scalar_v6_preflight_plan_20261004_v1.json')
    inventory_path = Path('metadata/baliphy_scalar_v6_preflight_launches_20261004_v1.json')
    plan = json.loads(plan_path.read_text()); verify(plan['pins'])
    inventory = json.loads(inventory_path.read_text())
    assert inventory['source_plan_sha256'] == sha(plan_path)
    assert inventory['launcher_sha256'] == sha('scripts/launch_baliphy_scalar_v6_preflight_v1.py')
    limits = {'cpu.max': '200000 100000', 'memory.max': str(32 * 2**30), 'memory.swap.max': '0'}
    handles = [observe(path, limits) for path in inventory['launches']]
    root = Path(plan['output'])
    rows = [json.loads(path.read_text()) for path in (root / 'chains').glob('*.json')]
    assert len({r['chain_id'] for r in rows}) == len(rows) <= 1620
    assert all(r['scientific_eligibility'] is False and r.get('posterior_sampling_launched', False) is False for r in rows)
    completed = None
    if Path(plan['completion']).exists():
        completed = json.loads(Path(plan['completion']).read_text())
        assert sha(completed['full_hash_archive']) == completed['full_hash_archive_sha256']
    pins = dict(plan['pins'])
    for path in [Path(__file__), plan_path, inventory_path,
                 Path('scripts/launch_baliphy_scalar_v6_preflight_v1.py')]: bind(pins, path)
    for path in inventory['launches']: bind(pins, path)
    result = dict(status='verified_original_scalar_v6_full_startup_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), source_plan=str(plan_path),
        source_plan_sha256=sha(plan_path), frozen_pins_checked=len(plan['pins']),
        original_handles=handles, unclosed_startup_checkpoints=len(rows),
        startup_status_counts=dict(Counter(r['status'] for r in rows)), expected_startup_roles=1620,
        producer_receipt_present=(root / 'receipt.json').exists(),
        serialized_readback_present=(root / 'readback.json').exists(), accounting_closure=completed,
        source_hashes=pins, scientific_eligibility=False, posterior_sampling_launched=False,
        posterior_qualified=False, gpu=False, all_eight_aims_incomplete=True,
        scope='Original exact controller/native identities, original journals and actual cgroup/native '
              'limits freshly observed. Frozen pins checked. Checkpoint counts remain unclosed '
              'progress; no full output replay, precise final native peaks, MCMC/posterior acceptance '
              'or biological inference. No original controller or failed attempt restarted.')
    with a.output.open('x') as f: json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes', 'original_handles', 'accounting_closure']}, indent=2))


if __name__ == '__main__':
    main()
