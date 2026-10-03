#!/usr/bin/env python3
"""Full 1,620-role independent native replay after original sampler closure."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil

import numpy as np

from ancestral_chain_attempt import sha, write_json
from independent_short_sampler_outputs_v2 import replay, summarize
from run_baliphy_reference_preflight import verify


PRODUCER_STATUS = 'complete_full_independent_short_sampler_native_replay_v2_pending_readback'
READER_STATUS = 'passed_full_independent_short_sampler_native_serialized_readback_v2'
COMPLETED_STATUS = 'complete_verified_full_independent_short_sampler_native_replay_v2'


def load(plan, path):
    verify(plan); bindings = dict(plan['pins']); bindings[str(path)] = sha(path)
    native_path = Path(plan['sampler_plan']); native = json.loads(native_path.read_text()); verify(native)
    completed_path = Path(native['completion']); completed = json.loads(completed_path.read_text())
    assert completed['status'] == 'complete_verified_full_reference_short_sampler_qualification'
    assert completed['full_chains'] == 1620 and completed['full_quartets'] == 405
    assert completed['checked_sampler_attempts'] + completed['unsuccessful_sampler_attempts'] == 1620
    assert completed['effective_inputs'] == 135 and completed['original_configuration_aliases'] == 324
    assert completed['posterior_qualified'] is completed['scientific_eligibility'] is False
    archive_path = Path(completed['full_hash_archive'])
    assert sha(archive_path) == completed['full_hash_archive_sha256']
    archive = json.loads(archive_path.read_text()); assert len(archive['services']) == 2
    for name, value in archive['source_hashes'].items():
        assert name not in bindings or bindings[name] == value; bindings[name] = value
    for p in [native_path, completed_path, archive_path]: bindings[str(p)] = sha(p)
    disposition_path = Path(native['output']) / 'dispositions.json'
    assert sha(disposition_path) == bindings[str(disposition_path)]
    dispositions = json.loads(disposition_path.read_text())
    jobs = json.loads(Path(native['jobs']).read_text())
    assert len(jobs) == len(dispositions) == 1620
    assert len({j['chain']['chain_id'] for j in jobs}) == len({r['chain_id'] for r in dispositions}) == 1620
    assert {j['chain']['chain_id'] for j in jobs} == {r['chain_id'] for r in dispositions}
    assert len({j['chain']['seed'] for j in jobs}) == 1620
    assert {j['chain']['seed'] for j in jobs}.isdisjoint({j['source_seed'] for j in jobs})
    verify({'pins': bindings})
    return jobs, {r['chain_id']: r for r in dispositions}, bindings


def check_array(path, result, states, free):
    with np.load(path, allow_pickle=False) as arrays:
        assert set(arrays.files) == {'states', 'iterations', 'unanchored_residue_counts'}
        assert arrays['states'].dtype == np.uint8 and np.array_equal(arrays['states'], states)
        assert arrays['unanchored_residue_counts'].dtype == np.int32
        assert np.array_equal(arrays['unanchored_residue_counts'], free)
        assert arrays['iterations'].dtype == np.int32 and arrays['iterations'].tolist() == [0, 10, 20]
    assert sha(path) == result['projection_array_sha256']


def run(path, reader=False, enforce_caps=True):
    plan = json.loads(path.read_text()); root = Path(plan['output'])
    if enforce_caps:
        group = next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
        cg = Path('/sys/fs/cgroup') / group.lstrip('/')
        limits = {k: (cg / k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
        assert limits == {'cpu.max': '200000 100000', 'memory.max': str(16 * 2**30), 'memory.swap.max': '0'}
        assert all(os.environ.get(k) == '1' for k in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'])
    jobs, dispositions, bindings = load(plan, path); initial = dict(bindings)
    root.mkdir(exist_ok=True); (root / 'chains').mkdir(exist_ok=True)
    lock = (root / ('readback.lock' if reader else 'stage.lock')).open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    final = root / ('readback.json' if reader else 'receipt.json')
    assert not final.exists(), 'Completed native replay refuses restart'
    marker = root / 'stage_plan.json'
    if reader or marker.exists(): assert json.loads(marker.read_text()) == plan
    else: write_json(marker, plan)
    producer = None
    if reader:
        producer_path = root / 'receipt.json'; producer = json.loads(producer_path.read_text())
        assert producer['status'] == PRODUCER_STATUS and producer['plan_sha256'] == sha(path)
        assert producer['source_hashes'] == initial
        for name, value in producer['artifacts'].items(): assert sha(root / name) == value
    rows = []; artifacts = {'stage_plan.json': sha(marker)}; expected_paths = set()
    for job in sorted(jobs, key=lambda j: j['chain']['chain_id']):
        assert shutil.disk_usage(root).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
        cid = job['chain']['chain_id']; checkpoint = root / 'chains' / (cid + '.json')
        array_path = checkpoint.with_suffix('.npz'); temporary = checkpoint.with_suffix('.npz.partial')
        existing = checkpoint.exists()
        assert not reader or existing, 'Missing whole-role checkpoint'
        assert not temporary.exists(), 'Incomplete array requires review'
        decoded = replay(job, dispositions[cid], plan['mapping'],
                         array_path=temporary if not reader and not existing and
                         dispositions[cid]['status'] == 'full_short_sampler_output_integrity_checked_not_posterior' else None)
        if isinstance(decoded, tuple):
            row, states, free = decoded
            if not reader and not existing:
                assert not array_path.exists(), 'Orphaned projection array requires review'
                os.replace(temporary, array_path)
            row.update(projection_array=str(array_path.relative_to(root)), projection_array_sha256=sha(array_path))
            check_array(array_path, row, states, free)
            expected_paths.add(array_path); artifacts[str(array_path.relative_to(root))] = sha(array_path)
        else:
            row = decoded
            assert not array_path.exists(), 'Unresolved role cannot acquire a successful array'
        if existing: assert json.loads(checkpoint.read_text()) == row, 'Changed serialized role'
        else: write_json(checkpoint, row)
        expected_paths.add(checkpoint); artifacts[str(checkpoint.relative_to(root))] = sha(checkpoint); rows.append(row)
        if len(rows) % 50 == 0: print('independent_short_native_replay', 'reader' if reader else 'producer', len(rows), '/1620', flush=True)
    assert {p for p in (root / 'chains').iterdir() if p.is_file()} == expected_paths, 'Foreign role artifact'
    summary = summarize(rows, jobs); export = root / 'dispositions.json'
    if reader:
        assert json.loads(export.read_text()) == rows
        assert all(producer[key] == value for key, value in summary.items())
    else: write_json(export, rows)
    artifacts['dispositions.json'] = sha(export)
    if reader: assert artifacts == producer['artifacts']
    verify({'pins': bindings})
    result = dict(status=READER_STATUS if reader else PRODUCER_STATUS, plan_sha256=sha(path), **summary,
        source_hashes=bindings, scientific_eligibility=False, third_independent_decoder=False, scope=plan['scope'])
    if reader:
        result['producer_receipt_sha256'] = sha(root / 'receipt.json')
        result['source_hashes'] = dict(bindings, **{str(root / 'receipt.json'): sha(root / 'receipt.json')},
                                     **{str(root / name): value for name, value in artifacts.items()})
    else: result['artifacts'] = artifacts
    write_json(final, result); print(json.dumps(summary), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True); parser.add_argument('--reader', action='store_true')
    args = parser.parse_args(); run(args.plan, args.reader)
