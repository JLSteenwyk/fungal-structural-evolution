#!/usr/bin/env python3
"""Run every joint-logger role after all four full-grid prerequisite closures.

This computational qualification is not a posterior-convergence check.
No role is dropped or automatically retried after a native failure.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
import json
import os
from pathlib import Path
import shutil
import time

import psutil

from ancestral_chain_attempt import live_group, run_attempt, sha, write_json
from baliphy_joint_sampler_qualification_v2 import inspect, summarize, SUMMARY_FIELDS
from baliphy_joint_sampler_gates import prerequisites
from reference_sampler_memory_budget import MemoryBudget
from run_baliphy_reference_preflight import verify


PRODUCER_STATUS = 'complete_full_joint_short_sampler_dispositions_pending_readback_v2'
READER_STATUS = 'passed_full_joint_short_sampler_serialized_readback_v2'
COMPLETED_STATUS = 'complete_verified_full_joint_short_sampler_qualification_v2'


def startup_gate(plan):
    return prerequisites(plan)


def verify_ledger(events, jobs, capacity, maximum_workers):
    """Independently reconstruct every admission/release and memory total."""
    expected = {j['chain']['chain_id']: j['memory_reservation_bytes'] for j in jobs}
    assert len(expected) == len(jobs)
    active = {}; acquired = set(); released = set(); total = peak = maximum_active = 0
    for i, event in enumerate(events):
        assert event['sequence'] == i
        identifier = event['identifier']; amount = event['amount']
        assert amount == expected[identifier]
        if event['action'] == 'acquire':
            assert identifier not in acquired
            acquired.add(identifier); active[identifier] = amount; total += amount
        else:
            assert event['action'] == 'release' and identifier not in released
            assert active.pop(identifier) == amount
            released.add(identifier); total -= amount
        assert 0 <= total <= capacity and len(active) <= maximum_workers
        assert event['reserved_total'] == total and event['active_roles'] == len(active)
        peak = max(peak, total); maximum_active = max(maximum_active, len(active))
    assert acquired == released == set(expected) and not active and total == 0
    return dict(roles_reserved=len(acquired), ledger_events=len(events), peak_reserved_bytes=peak,
                maximum_simultaneous_roles=maximum_active, reservation_capacity_bytes=capacity,
                scope='Declared address-space reservations held through native execution and output audit; not measured native peak RSS.')


def runtime_caps(plan, reader=False):
    resources = plan['resources']; cpus = 2 if reader else resources['cpus']
    memory = 32 if reader else resources['memory_gib']
    group = next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
    root = Path('/sys/fs/cgroup') / group.lstrip('/')
    limits = {k: (root / k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert limits == {'cpu.max': str(cpus * 100000) + ' 100000',
                      'memory.max': str(memory * 2**30), 'memory.swap.max': '0'}
    assert all(os.environ.get(k) == '1' for k in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'])
    return limits


def execute_job(job, root, digest, mapping, budget, minimum_free_disk_gib):
    """A native exit failure is an outcome; controller/integrity errors abort admission."""
    cid = job['chain']['chain_id']; folder = root / 'attempts' / cid
    with budget.reserve(cid, job['memory_reservation_bytes']):
        try:
            assert shutil.disk_usage(root).free >= minimum_free_disk_gib * 2**30
            previous = sorted(folder.glob('attempt-[0-9]*'))
            if previous:
                assert len(previous) == 1, 'Never silently choose or retry qualification attempts'
                receipt = previous[0] / 'receipt.json'
                assert receipt.is_file(), 'Incomplete native attempt requires manual review'
                identity = json.loads((previous[0] / 'process.json').read_text())
                assert not live_group(identity['pgid']), 'Original native group is still active'
            else:
                receipt = run_attempt(folder, job['config'])
            row = inspect(job, receipt, digest, mapping, root / "frames" / cid)
            cp = root / 'chains' / (cid + '.json')
            if cp.exists(): assert json.loads(cp.read_text()) == row
            else: write_json(cp, row)
            return row
        except BaseException as error:
            # A failed output check must not release memory and admit another
            # native job if the exceptional attempt still has live descendants.
            budget.abort(cid + ': ' + type(error).__name__ + ': ' + str(error))
            raise


def run(path, reader=False):
    plan = json.loads(path.read_text()); verify(plan); digest = sha(path)
    limits = runtime_caps(plan, reader); bindings = dict(plan['pins'])
    bindings[str(path)] = digest
    for p, h in startup_gate(plan).items():
        assert p not in bindings or bindings[p] == h; bindings[p] = h
    jobs = json.loads(Path(plan['jobs']).read_text()); assert len(jobs) == 1620
    root = Path(plan['output']); root.mkdir(exist_ok=True)
    lock = (root / ('readback.lock' if reader else 'stage.lock')).open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    final = root / ('readback.json' if reader else 'receipt.json')
    assert not final.exists(), 'Completed qualification cannot restart'
    (root / 'chains').mkdir(exist_ok=True); journal = root / 'memory_reservations.jsonl'
    started = time.monotonic(); rows = []
    if reader:
        producer_path = root / 'receipt.json'; producer = json.loads(producer_path.read_text())
        assert producer['status'] == PRODUCER_STATUS and producer['plan_sha256'] == digest
        verify({'pins': producer['source_hashes']})
        for job in jobs:
            cp = root / 'chains' / (job['chain']['chain_id'] + '.json')
            saved = json.loads(cp.read_text())
            reconstructed = inspect(job, Path(saved['native_receipt']), digest, plan['mapping'], root / 'frames' / job['chain']['chain_id'], allow_export_creation=False)
            assert saved == reconstructed; rows.append(reconstructed)
            if len(rows) % 50 == 0: print('joint_sampler_readback', len(rows), '/1620', flush=True)
        rows.sort(key=lambda x: x['chain_id'])
        assert json.loads((root / 'dispositions.json').read_text()) == rows
        assert json.loads((root / 'stage_plan.json').read_text()) == plan
        for name, h in producer['artifacts'].items(): assert sha(root / name) == h
    else:
        assert psutil.virtual_memory().available >= plan['resources']['memory_gib'] * 2**30
        assert not journal.exists(), 'Partial controller run requires manual review; never silently resume reservations'
        marker = root / 'stage_plan.json'
        if marker.exists(): assert json.loads(marker.read_text()) == plan
        else: write_json(marker, plan)
        with journal.open('x') as handle:
            def persist(event):
                handle.write(json.dumps(event, allow_nan=False) + '\n'); handle.flush()
            budget = MemoryBudget(plan['resources']['reservation_capacity_gib'] * 2**30, persist)
            with ThreadPoolExecutor(max_workers=plan['resources']['workers']) as pool:
                for future in as_completed([pool.submit(execute_job, j, root, digest, plan['mapping'], budget,
                    plan['resources']['minimum_free_disk_gib']) for j in jobs]):
                    row = future.result(); rows.append(row)
                    print('joint_sampler_qualifications', len(rows), '/1620', row['chain_id'], row['status'], flush=True)
            assert budget.aborted is None
        rows.sort(key=lambda x: x['chain_id']); write_json(root / 'dispositions.json', rows)
    events = [json.loads(x) for x in journal.read_text().splitlines()]
    ledger = verify_ledger(events, jobs, plan['resources']['reservation_capacity_gib'] * 2**30, plan['resources']['workers'])
    summary = summarize(rows)
    for r in rows:
        p = Path(r['native_receipt']); bindings[str(p)] = sha(p)
        bindings[str(p.parent.parent / 'configuration.json')] = sha(p.parent.parent / 'configuration.json')
        for name, h in json.loads(p.read_text())['artifacts'].items(): bindings[str(p.parent / name)] = h
        for frame in r['joint_frames']:
            assert sha(frame['projection_array']) == frame['projection_array_sha256']
            bindings[frame['projection_array']] = frame['projection_array_sha256']
    expected_exports = {Path(frame['projection_array']) for row in rows for frame in row['joint_frames']}
    assert set((root / 'frames').glob('*/*')) == expected_exports, 'Foreign joint frame artifact'
    verify(dict(pins=bindings)); assert sha(path) == digest
    result = dict(status=READER_STATUS if reader else PRODUCER_STATUS, plan_sha256=digest,
        **summary, reservation_audit=ledger, actual_cgroup_limits=limits, source_hashes=bindings,
        elapsed_seconds=time.monotonic() - started, scientific_eligibility=False, scope=plan['scope'])
    if reader:
        assert all(producer[k] == result[k] for k in SUMMARY_FIELDS + ['reservation_audit'])
        result['producer_receipt_sha256'] = sha(producer_path)
    else:
        result['artifacts'] = {str(p.relative_to(root)): sha(p) for p in [root / 'stage_plan.json',
            root / 'dispositions.json', journal, *sorted((root / 'chains').glob('*.json'))]}
    with final.open('x') as handle: handle.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(summary), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--reader', action='store_true'); args = parser.parse_args(); run(args.plan, args.reader)
