#!/usr/bin/env python3
"""Recalculate the full closed recovery grid with separate direct-lag estimators."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import time

from independent_ancestral_scalar_diagnostics import diagnose, compare
from independent_baliphy_scalar_sources import load, trace_input, summarize, numeric_status, CUTS, KINDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def atomic(path, value):
    temporary = path.with_suffix(path.suffix + '.partial')
    with temporary.open('w') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n')
        f.flush(); os.fsync(f.fileno())
    os.replace(temporary, path)


def group_result(source, bindings, group, info, plan):
    rows = []; cache = {}
    if info['status'] != 'unresolved_failed_native_chain_retained':
        for kind in KINDS:
            for cutoff in CUTS:
                variables, values, report, output = trace_input(source, bindings, group, info, kind, cutoff, cache)
                for i, variable in enumerate(variables):
                    independent = diagnose(values[:, :, i]); original = report['variables'][variable]
                    errors, unresolved = compare(original, independent, values=values[:, :, i],
                        **plan['comparison_tolerances'])
                    row = dict(kind=kind, cutoff=cutoff, variable=variable,
                        original_report=output['path'], original_report_sha256=output['sha256'],
                        original=original, independent=independent, absolute_errors=errors,
                        unresolved_numeric_metrics=unresolved, scientific_eligibility=False)
                    row['numeric_status'] = numeric_status(row); rows.append(row)
    return dict(status=('full_scalar_length_numerical_comparison_not_posterior_qualification'
        if rows else 'unresolved_failed_native_chain_retained'), chain_ids=info['chain_ids'],
        rows=rows, scientific_eligibility=False)


def run(path, stop_after_groups=None):
    started = time.perf_counter(); plan = json.loads(path.read_text())
    source, bindings = load(plan, path); root = Path(plan['output'])
    assert shutil.disk_usage(root.parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    root.mkdir(exist_ok=True)
    lock = (root / 'stage.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (root / 'receipt.json').exists(), 'Completed stage cannot restart'
    stage = dict(schema='independent-baliphy-scalar-direct-lag-v1', plan_sha256=sha(path),
        diagnostic_archive_sha256=source['complete']['full_hash_archive_sha256'])
    marker = root / 'stage_plan.json'
    if marker.exists():
        assert json.loads(marker.read_text()) == stage
    else:
        atomic(marker, stage)
    folder = root / 'groups'; folder.mkdir(exist_ok=True); results = {}; manifest = []
    for number, (group, info) in enumerate(sorted(source['groups'].items()), 1):
        cp = folder / (group + '.json')
        expected = group_result(source, bindings, group, info, plan)
        if cp.exists():
            saved = json.loads(cp.read_text()); assert saved == dict(stage=stage, group=group, result=expected)
        else:
            atomic(cp, dict(stage=stage, group=group, result=expected))
        results[group] = expected
        manifest.append(dict(group=group, path=str(cp.relative_to(root)), sha256=sha(cp)))
        if number % 25 == 0:
            print('independent_scalar_quartets', number, '/', len(source['groups']), flush=True)
        if stop_after_groups == number:
            raise InterruptedError('Software checkpoint contract')
    summary = summarize(results, source); atomic(root / 'group_manifest.json', manifest)
    verify(bindings)
    artifacts = {str(p.relative_to(root)): sha(p) for p in [marker, root / 'group_manifest.json', *[root / m['path'] for m in manifest]]}
    result = dict(status='complete_full_independent_baliphy_scalar_comparison_pending_readback',
        plan_sha256=sha(path), **summary, source_hashes=bindings, artifacts=artifacts,
        elapsed_seconds=time.perf_counter() - started, scientific_eligibility=False, scope=plan['scope'])
    atomic(root / 'receipt.json', result); print(json.dumps(summary), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args(); run(args.plan)
