#!/usr/bin/env python3
"""Rebuild every serialized raw-trace marginal comparison in the full grid."""
import argparse
import fcntl
import json
from pathlib import Path

from independent_baliphy_scalar_sources import load, summarize
from prepare_independent_baliphy_scalar_readback_v2 import group_result
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(path, output):
    plan = json.loads(path.read_text()); root = Path(plan['output'])
    lock = (root / 'readback.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (root / 'readback.json').exists(), 'Completed original reader cannot restart'
    source, bindings = load(plan, path); original = dict(bindings)
    rp = root / 'receipt.json'; receipt = json.loads(rp.read_text()); bind(bindings, rp)
    assert receipt['status'] == 'complete_full_independent_baliphy_scalar_comparison_pending_readback'
    assert receipt['plan_sha256'] == sha(path) and receipt['source_hashes'] == original
    assert receipt['scientific_eligibility'] is False
    for name, digest in receipt['artifacts'].items():
        bind(bindings, root / name, digest)
    verify(bindings)
    stage = dict(schema='independent-baliphy-scalar-direct-lag-v2', plan_sha256=sha(path),
        diagnostic_archive_sha256=source['complete']['full_hash_archive_sha256'])
    assert json.loads((root / 'stage_plan.json').read_text()) == stage
    manifest = json.loads((root / 'group_manifest.json').read_text())
    assert [m['group'] for m in manifest] == sorted(source['groups'])
    results = {}
    for group, entry in zip(sorted(source['groups']), manifest):
        cp = root / entry['path']; assert sha(cp) == entry['sha256']
        assert entry['path'] == 'groups/' + group + '.json'
        expected = group_result(source, bindings, group, source['groups'][group], plan)
        assert json.loads(cp.read_text()) == dict(stage=stage, group=group, result=expected)
        results[group] = expected
    summary = summarize(results, source)
    assert all(receipt[k] == v for k, v in summary.items())
    verify(bindings)
    result = dict(status='passed_full_independent_baliphy_scalar_serialized_comparison_readback',
        plan_sha256=sha(path), producer_receipt_sha256=sha(rp), **summary,
        source_hashes=bindings, scientific_eligibility=False,
        third_independent_metric_implementation=False, scope=plan['scope'])
    with Path(output).open('x') as f:
        f.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(summary), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True); args = parser.parse_args(); run(args.plan, args.output)
