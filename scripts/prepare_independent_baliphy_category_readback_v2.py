#!/usr/bin/env python3
"""Full recovered categorical numerical replay with immutable group checkpoints."""
import argparse
import fcntl
import json
from pathlib import Path
import shutil
import time

from independent_baliphy_category_sources_v2 import load, summarize
from independent_baliphy_category_replay_v2 import replay_group
from prepare_independent_baliphy_scalar_readback_v2 import atomic
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def stage_marker(path, source):
    return dict(schema='independent-baliphy-categorical-direct-lag-v2', plan_sha256=sha(path),
        diagnostic_archive_sha256=source['complete']['full_hash_archive_sha256'])


def artifacts(root, marker, manifest, entries, results):
    paths = [marker, manifest, *(root / e['path'] for e in entries)]
    paths.extend(root / c['serialized_comparison'] for r in results.values() for c in r['cutoffs'].values())
    assert len(paths) == len(set(paths))
    assert set((root / 'groups').glob('*.json')) == {root / e['path'] for e in entries}
    assert set((root / 'groups').rglob('*.jsonl.gz')) == {p for p in paths if p.name.endswith('.jsonl.gz')}
    return {str(p.relative_to(root)):sha(p) for p in paths}


def run(path, stop_after_groups=None):
    started = time.perf_counter(); plan = json.loads(path.read_text()); root = Path(plan['output'])
    assert not (root / 'receipt.json').exists(), 'Completed producer cannot restart'
    assert shutil.disk_usage(root.parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    root.mkdir(exist_ok=True)
    lock = (root / 'stage.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (root / 'receipt.json').exists(), 'Completed producer cannot restart'
    source, bindings = load(plan, path)
    stage = stage_marker(path, source); marker = root / 'stage_plan.json'
    if marker.exists(): assert json.loads(marker.read_text()) == stage
    else: atomic(marker, stage)
    folder = root / 'groups'; folder.mkdir(exist_ok=True); results = {}; entries = []
    for number, (group, info) in enumerate(sorted(source['groups'].items()), 1):
        cp = folder / (group + '.json')
        expected = replay_group(source, bindings, group, info, root, plan)
        checkpoint = dict(stage=stage, group=group, result=expected)
        if cp.exists(): assert json.loads(cp.read_text()) == checkpoint
        else: atomic(cp, checkpoint)
        results[group] = expected
        entries.append(dict(group=group, path=str(cp.relative_to(root)), sha256=sha(cp)))
        print('independent_categorical_quartets', number, '/', len(source['groups']),
            'elapsed_seconds', round(time.perf_counter() - started, 3), flush=True)
        if stop_after_groups == number: raise InterruptedError('Software checkpoint contract')
    summary = summarize(results, source); manifest = root / 'group_manifest.json'; atomic(manifest, entries)
    verify(bindings)
    exported = artifacts(root, marker, manifest, entries, results)
    result = dict(status='complete_full_independent_baliphy_categorical_comparison_pending_readback',
        plan_sha256=sha(path), **summary, source_hashes=bindings, artifacts=exported,
        elapsed_seconds=time.perf_counter() - started, scientific_eligibility=False, scope=plan['scope'])
    atomic(root / 'receipt.json', result); print(json.dumps(summary), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args(); run(args.plan)
