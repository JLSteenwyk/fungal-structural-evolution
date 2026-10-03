#!/usr/bin/env python3
"""Check selected immutable completed draws while the full producer continues."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import psutil

from ancestral_chain_attempt import sha
from matched_predictor_branch_inputs import verify
from matched_predictor_resampling import load, add_splits
from readback_matched_predictor_resampling import read_case
from record_project_runtime_checkpoint_v4 import fingerprint


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--scratch', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists() and not args.scratch.exists()
    args.scratch.mkdir(parents=True)
    path = Path('metadata/matched_predictor_resampling_plan_20261003_v1.json')
    plan = json.loads(path.read_text()); source, bindings = load(plan, path)
    add_splits(source, Path(plan['point_output'])); root = Path(plan['output'])
    inventory = json.loads(Path(plan['launch_inventory']).read_text())
    launch_path = Path(inventory['launches'][0]); launch = json.loads(launch_path.read_text())
    process = fingerprint(launch); assert process is not None
    completed = sorted((root / 'cases').glob('*/*/*.json')); assert completed
    # Check both ends of every available input/mode; at most twelve cases.
    groups = {}
    for item in completed: groups.setdefault((item.parent.parent.name, item.parent.name), []).append(item)
    selected = []
    for paths in groups.values():
        for item in [paths[0], paths[-1]]:
            if item not in selected and len(selected) < 12: selected.append(item)
    native_workers = []
    deadline = time.monotonic() + 10
    while not native_workers and time.monotonic() < deadline:
        assert fingerprint(launch) is not None
        for child in process.children(recursive=True):
            try:
                command = child.cmdline()
                if not command or command[0] != plan['executable']: continue
                limits = dict(address_space=child.rlimit(psutil.RLIMIT_AS), cpu=child.rlimit(psutil.RLIMIT_CPU), file=child.rlimit(psutil.RLIMIT_FSIZE))
                assert limits == dict(address_space=(2 * 2**30, 2 * 2**30), cpu=(300, 300), file=(8 * 2**20, 8 * 2**20))
                native_workers.append(dict(pid=child.pid, created=child.create_time(), command=command, limits=limits))
            except psutil.NoSuchProcess: continue
        if not native_workers: time.sleep(.05)
    rows = []
    for item in selected:
        saved = json.loads(item.read_text()); key = saved['original_input_id']
        result = read_case(plan, source, key, saved['mode'], saved['replicate'], root, args.scratch)
        assert result == saved
        for field in ['archive', 'array']:
            bindings[saved[field]] = saved[field + '_sha256']
        bindings[str(item)] = sha(item); rows.append(saved)
    verify(bindings); assert not list(args.scratch.iterdir())
    result = dict(status='passed_partial_available_native_resampling_readback', checked_utc=datetime.now(timezone.utc).isoformat(),
        immutable_completed_cases_at_snapshot=len(completed), checked_cases=len(rows), checked_native_roles=len(rows) * 7,
        checked_modes=sorted({row['mode'] for row in rows}), checked_original_inputs=len({row['original_input_id'] for row in rows}),
        checked_status_counts=dict(sum((Counter(row['native_status_counts']) for row in rows), Counter())),
        original_producer_pid=process.pid, current_sampled_native_worker_limits=native_workers,
        source_hashes=bindings, scientific_eligibility=False,
        scope='Read-only independent archive/tree/report/array checks of at most twelve completed production cases, explicitly partial. Full53200-case producer, reader and closure continue unchanged. No native runs, retries, completed-output writes or scientific confidence intervals.')
    result['source_hashes'][str(Path(__file__))] = sha(__file__)
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ['source_hashes', 'current_sampled_native_worker_limits']}, indent=2))
    print('sampled_native_workers', len(native_workers))


if __name__ == '__main__': main()
