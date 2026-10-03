#!/usr/bin/env python3
"""Replay every native startup, recognizing only the native timing footer.

Original programs, plans, native attempts and v1 dispositions stay immutable.
This reads exactly the same full grid; it never launches BAli-Phy.
"""
import argparse
from datetime import datetime
import fcntl
import json
import math
from pathlib import Path
import re

from ancestral_chain_attempt import sha, write_json
from baliphy_reference_initialization import alignment_records, transform, validate_initial_json
from run_baliphy_reference_preflight import inspect, summarize, tree_check, verify


def parse_initial_model(text):
    """One JSON model plus whitespace or the exact five-line native Work footer."""
    stripped = text.lstrip(); value, end = json.JSONDecoder().raw_decode(stripped)
    assert isinstance(value, dict) and value['iter'] == 0
    rest = stripped[end:].strip()
    if not rest: return value, None
    lines = rest.splitlines()
    assert len(lines) == 5 and lines[0] == 'Work:', 'Unexpected native stdout after initial model'
    dates = []
    for line, name in zip(lines[1:3], ['start', 'end']):
        match = re.fullmatch(r'\s*' + name + r':\s*(.+)', line); assert match
        dates.append(datetime.strptime(match[1], '%a %b %d %H:%M:%S %Y'))
    assert dates[1] >= dates[0]
    seconds = []
    for line, name in zip(lines[3:], ['elapsed', 'CPU']):
        match = re.fullmatch(r'\s*total \(' + name + r'\) time:\s*'
                            r'[0-9]+[dhms](?:\s+[0-9]+[dhms])*\s+\(([0-9]+(?:\.[0-9]+)?)s\)', line)
        assert match, 'Unexpected native timing format'
        number = float(match[1]); assert math.isfinite(number) and number >= 0; seconds.append(number)
    return value, dict(start_local=dates[0].isoformat(), end_local=dates[1].isoformat(),
                       reported_elapsed_seconds=seconds[0], reported_cpu_seconds=seconds[1],
                       scope='Native self-reported Work footer; not whole-cgroup peak memory or convergence evidence')


def inspect_v2(job, receipt, native_digest, stage_digest):
    original = inspect(job, receipt, native_digest)
    row = {k: v for k, v in original.items() if k not in ['status', 'audit', 'error_type', 'error', 'plan_sha256', 'posterior_sampling_launched']}
    row.update(original_startup_disposition=original['status'], native_plan_sha256=native_digest,
               plan_sha256=stage_digest)
    if original['exit_code'] != 0: return dict(**row, status='unsuccessful_native_startup_retained')
    try:
        chain = job['chain']
        assert Path(job['program']).read_text() == transform(Path(chain['program']).read_text())
        value, timing = parse_initial_model((receipt.parent / 'stdout.log').read_text())
        tips, nodes = tree_check(chain['tree'], receipt.parent / 'runtime-tree.nwk')
        audit = validate_initial_json(value, alignment_records(chain['alignment']), tips)
        assert audit['nodes'] == nodes and audit['tips'] == chain['proteins']
        return dict(**row, status='reference_startup_homology_density_and_representation_checked',
                    audit=audit, stdout_timing_summary=timing, posterior_sampling_launched=False)
    except (AssertionError, KeyError, ValueError) as error:
        return dict(**row, status='invalid_native_startup_retained',
                    error_type=type(error).__name__, error=str(error))


def reconstruct(path):
    plan = json.loads(path.read_text()); verify(plan); bindings = dict(plan['pins']); bindings[str(path)] = sha(path)
    native_path = Path(plan['native_plan']); native = json.loads(native_path.read_text()); verify(native)
    closed_path = Path(native['completion']); closed = json.loads(closed_path.read_text())
    assert closed['status'] == 'complete_verified_full_reference_alignment_startup'
    archive = Path(closed['full_hash_archive']); assert sha(archive) == closed['full_hash_archive_sha256']
    proof = json.loads(archive.read_text()); assert len(proof['services']) == 2
    for p, h in proof['source_hashes'].items():
        assert p not in bindings or bindings[p] == h; bindings[p] = h
    bindings[str(closed_path)] = sha(closed_path); bindings[str(archive)] = sha(archive)
    native_root = Path(native['output']); native_receipt = json.loads((native_root / 'receipt.json').read_text())
    assert native_receipt['plan_sha256'] == sha(native_path)
    rows = []; jobs = json.loads(Path(native['jobs']).read_text()); assert len(jobs) == 1620
    for number, job in enumerate(sorted(jobs, key=lambda x: x['chain']['chain_id']), 1):
        cp = native_root / 'chains' / (job['chain']['chain_id'] + '.json'); saved = json.loads(cp.read_text())
        receipt = Path(saved['attempt_receipt']); original = inspect(job, receipt, sha(native_path))
        assert saved == original
        rows.append(inspect_v2(job, receipt, sha(native_path), sha(path)))
        if number % 50 == 0: print('reference_startup_footer_replay', number, '/1620', flush=True)
    verify(dict(pins=bindings))
    summary = summarize(rows)
    summary['native_timing_footers_checked'] = sum(r.get('stdout_timing_summary') is not None for r in rows)
    summary['v1_parser_dispositions_reclassified'] = sum(r['status'] != r['original_startup_disposition'] for r in rows)
    return plan, rows, summary, bindings


def run(path, reader=False):
    plan = json.loads(path.read_text()); root = Path(plan['output']); root.mkdir(exist_ok=True)
    lock = (root / ('readback.lock' if reader else 'stage.lock')).open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    final = root / ('readback.json' if reader else 'receipt.json')
    assert not final.exists(), 'Completed footer replay cannot restart'
    plan, rows, summary, bindings = reconstruct(path)
    exported = root / 'dispositions.json'
    if reader:
        producer_path = root / 'receipt.json'; producer = json.loads(producer_path.read_text())
        assert producer['status'] == 'complete_full_reference_startup_footer_replay_pending_readback'
        assert producer['plan_sha256'] == sha(path)
        assert json.loads(exported.read_text()) == rows
        assert all(producer[k] == v for k, v in summary.items())
        for name, h in producer['artifacts'].items(): assert sha(root / name) == h
        result = dict(status='passed_full_reference_startup_footer_serialized_readback',
            producer_receipt_sha256=sha(producer_path), plan_sha256=sha(path), **summary,
            source_hashes=bindings, scientific_eligibility=False, scope=plan['scope'])
    else:
        write_json(exported, rows)
        result = dict(status='complete_full_reference_startup_footer_replay_pending_readback',
            plan_sha256=sha(path), **summary, source_hashes=bindings,
            artifacts={'dispositions.json': sha(exported)}, scientific_eligibility=False, scope=plan['scope'])
    with final.open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(summary), flush=True); return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--reader', action='store_true'); a = p.parse_args(); run(a.plan, a.reader)
