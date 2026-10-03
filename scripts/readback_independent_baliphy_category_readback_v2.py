#!/usr/bin/env python3
"""Rebuild every complete categorical serialized comparison from source arrays."""
import argparse
import fcntl
import json
from pathlib import Path

from independent_baliphy_category_sources_v2 import load, summarize
from independent_baliphy_category_replay_v2 import replay_group
from prepare_independent_baliphy_category_readback_v2 import stage_marker, artifacts
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(path, output):
    plan = json.loads(path.read_text()); root = Path(plan['output'])
    lock = (root / 'readback.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (root / 'readback.json').exists(), 'Completed original reader cannot restart'
    source, bindings = load(plan, path); original_bindings = dict(bindings)
    rp = root / 'receipt.json'; receipt = json.loads(rp.read_text()); bind(bindings, rp)
    assert receipt['status'] == 'complete_full_independent_baliphy_categorical_comparison_pending_readback'
    assert receipt['plan_sha256'] == sha(path) and receipt['source_hashes'] == original_bindings
    assert receipt['scientific_eligibility'] is False
    for name, digest in receipt['artifacts'].items(): bind(bindings, root / name, digest)
    verify(bindings)
    stage = stage_marker(path, source); marker = root / 'stage_plan.json'
    assert json.loads(marker.read_text()) == stage
    manifest = root / 'group_manifest.json'; entries = json.loads(manifest.read_text())
    assert [e['group'] for e in entries] == sorted(source['groups'])
    results = {}
    for number, (group, entry) in enumerate(zip(sorted(source['groups']), entries), 1):
        assert entry['path'] == 'groups/' + group + '.json'
        cp = root / entry['path']; assert sha(cp) == entry['sha256']
        expected = replay_group(source, bindings, group, source['groups'][group], root, plan, readback=True)
        assert json.loads(cp.read_text()) == dict(stage=stage, group=group, result=expected)
        results[group] = expected
        print('readback_categorical_quartets', number, '/', len(source['groups']), flush=True)
    summary = summarize(results, source)
    assert all(receipt[k] == v for k,v in summary.items())
    assert artifacts(root, marker, manifest, entries, results) == receipt['artifacts']
    verify(bindings)
    result = dict(status='passed_full_independent_baliphy_categorical_serialized_comparison_readback',
        plan_sha256=sha(path), producer_receipt_sha256=sha(rp), **summary, source_hashes=bindings,
        scientific_eligibility=False, third_independent_metric_implementation=False, scope=plan['scope'])
    with Path(output).open('x') as f: f.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(summary), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True); args = parser.parse_args(); run(args.plan, args.output)
