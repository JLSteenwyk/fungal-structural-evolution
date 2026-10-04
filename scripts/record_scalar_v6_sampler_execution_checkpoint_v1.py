#!/usr/bin/env python3
"""Observe six original sampler/telemetry controllers and unclosed progress."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); assert not args.output.exists()
    plan_path = Path('metadata/baliphy_scalar_v6_sampler_execution_plan_20261004_v1.json')
    plan = json.loads(plan_path.read_text()); verify(plan['pins'])
    inventory_path = Path(plan['launch_inventory'])
    inventory = json.loads(inventory_path.read_text())
    assert inventory['source_plan_sha256'] == sha(plan_path)
    assert inventory['expected_native_roles'] == 1620
    assert inventory['iterations_per_role'] == 20 and not inventory['posterior_qualified']
    launches = inventory['sampler_launches'] + inventory['observer_launches']
    assert len(launches) == len(set(launches)) == 6
    expected = [(16, 200), (2, 32), (2, 32), (1, 2), (2, 4), (2, 4)]
    handles = [observe(path, {'cpu.max': str(cpu*100000)+' 100000',
        'memory.max': str(memory*2**30), 'memory.swap.max': '0'})
        for path, (cpu, memory) in zip(launches, expected)]
    prerequisites = [observe(path) for path in plan['dependencies']]
    root = Path(plan['output']); rows = []; transitioning = []
    for path in (root/'chains').glob('*.json'):
        try:
            rows.append(json.loads(path.read_text()))
        except json.JSONDecodeError:
            transitioning.append(str(path))
    assert len({row['chain_id'] for row in rows}) == len(rows) <= 1620
    assert all(row['scientific_eligibility'] is False for row in rows)
    receipts = {}
    for label, directory in [('sampler', root), ('observer', Path(plan['observer_output']))]:
        receipts[label] = {name: (directory/name).is_file() for name in ['receipt.json', 'readback.json']}
    closures = {}
    for key in ['completion', 'observer_completion']:
        path = Path(plan[key]); item = dict(path=str(path), present=path.exists())
        if path.exists():
            receipt = json.loads(path.read_text())
            assert sha(receipt['full_hash_archive']) == receipt['full_hash_archive_sha256']
            item.update(status=receipt['status'], sha256=sha(path))
        closures[key] = item
    pins = dict(plan['pins'])
    for path in [Path(__file__), plan_path, inventory_path, *map(Path, launches)]:
        bind(pins, path)
    result = dict(status='verified_original_scalar_v6_sampler_and_observer_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), plan=str(plan_path), plan_sha256=sha(plan_path),
        frozen_pins_checked=len(plan['pins']), original_handles=handles,
        original_prerequisite_handles=prerequisites, unclosed_disposition_rows=len(rows),
        disposition_status_counts=dict(Counter(row['status'] for row in rows)),
        transitioning_checkpoint_files=transitioning,
        original_attempt_identity_files=len(list((root/'attempts').glob('*/attempt-*/process.json'))),
        receipt_states=receipts, completion_gates=closures, source_hashes=pins,
        posterior_qualified=False, scientific_eligibility=False, biological_fits=0,
        gpu=False, all_eight_aims_incomplete=True,
        scope='Six exact original controller/native identities and cgroup limits, original prerequisites '
              'and unclosed checkpoint counts observed. Partially serialized files are explicit. '
              'No full replay, guaranteed final native peaks, posterior or biological acceptance. '
              'No failed attempt or original controller restarted.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in
        ['source_hashes', 'original_handles', 'original_prerequisite_handles']}, indent=2))


if __name__ == '__main__':
    main()
