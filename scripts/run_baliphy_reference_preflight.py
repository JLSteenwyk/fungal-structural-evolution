#!/usr/bin/env python3
"""Native --test for the entire frozen input/prior/seed grid, without MCMC."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

from Bio import Phylo

from ancestral_chain_attempt import run_attempt, sha, write_json
from baliphy_reference_initialization import alignment_records, transform, validate_initial_json


SUMMARY_FIELDS = ['full_chains', 'full_quartets', 'effective_inputs', 'original_configuration_aliases',
                  'validated_startups', 'unsuccessful_startups', 'complete_startup_quartets',
                  'unresolved_startup_quartets', 'status_counts', 'posterior_sampling_launched']


def verify(plan):
    for path, h in plan['pins'].items(): assert sha(path) == h, path


def tree_check(source_path, runtime_path):
    source, runtime = [Phylo.read(p, 'newick') for p in [source_path, runtime_path]]
    def nodes(tree):
        result = {}
        for n in tree.find_clades():
            key = tuple(sorted(x.name for x in n.get_terminals()))
            assert key not in result, 'Duplicate rooted descendant clade'
            result[key] = n
        return result
    first, second = nodes(source), nodes(runtime); assert first.keys() == second.keys()
    for key, node in first.items():
        if node is not source.root:
            assert abs(node.branch_length - second[key].branch_length) < 1e-10
    return [n.name for n in runtime.get_terminals()], len(second)


def inspect(job, receipt, plan_digest):
    attempt = json.loads(receipt.read_text())
    encoded = json.dumps(job['config'], sort_keys=True, separators=(',', ':'), allow_nan=False)
    assert attempt['configuration_sha256'] == hashlib.sha256(encoded.encode()).hexdigest()
    assert json.loads((receipt.parent.parent / 'configuration.json').read_text()) == job['config']
    assert json.loads((receipt.parent / 'command.json').read_text()) == job['config']['command']
    assert json.loads((receipt.parent / 'process.json').read_text())['command'] == job['config']['command']
    for name, h in attempt['artifacts'].items(): assert sha(receipt.parent / name) == h, name
    chain = job['chain']
    row = dict(chain_id=chain['chain_id'], effective_input_group=chain['effective_input_group'],
        model_input_identity=chain['effective_input_group'] + '-' + chain['prior_label'],
        chain_role=chain['chain'], prior_label=chain['prior_label'], family=chain['family'],
        original_configuration_ids=chain['original_configuration_ids'], source_seed=chain['seed'],
        fresh_seed=job['fresh_seed'], attempt_receipt=str(receipt), attempt_receipt_sha256=sha(receipt),
        exit_code=attempt['exit_code'], native_status=attempt['status'],
        elapsed_seconds=attempt['elapsed_seconds'], plan_sha256=plan_digest,
        scientific_eligibility=False)
    if attempt['exit_code'] != 0:
        return dict(**row, status='unsuccessful_native_startup_retained')
    try:
        assert Path(job['program']).read_text() == transform(Path(chain['program']).read_text())
        value = json.loads((receipt.parent / 'stdout.log').read_text())
        assert value['iter'] == 0
        tips, nodes = tree_check(chain['tree'], receipt.parent / 'runtime-tree.nwk')
        audit = validate_initial_json(value, alignment_records(chain['alignment']), tips)
        assert audit['nodes'] == nodes and audit['tips'] == chain['proteins']
        return dict(**row, status='reference_startup_homology_density_and_representation_checked',
                    audit=audit, posterior_sampling_launched=False)
    except (AssertionError, KeyError, ValueError) as error:
        return dict(**row, status='invalid_native_startup_retained',
                    error_type=type(error).__name__, error=str(error))


def summarize(rows):
    assert len(rows) == len({r['chain_id'] for r in rows}) == 1620
    groups = defaultdict(list)
    for r in rows: groups[r['model_input_identity']].append(r)
    assert len(groups) == 405
    assert all(len(v) == 4 and {r['chain_role'] for r in v} == {1, 2, 3, 4} for v in groups.values())
    valid = 'reference_startup_homology_density_and_representation_checked'
    complete = sum(all(r['status'] == valid for r in v) for v in groups.values())
    return dict(full_chains=len(rows), full_quartets=len(groups),
        effective_inputs=len({r['effective_input_group'] for r in rows}),
        original_configuration_aliases=len({x for r in rows for x in r['original_configuration_ids']}),
        validated_startups=sum(r['status'] == valid for r in rows),
        unsuccessful_startups=sum(r['status'] != valid for r in rows),
        complete_startup_quartets=complete, unresolved_startup_quartets=405 - complete,
        status_counts=dict(Counter(r['status'] for r in rows)), posterior_sampling_launched=False)


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args(); plan = json.loads(a.plan.read_text()); digest = sha(a.plan); verify(plan)
    group = next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
    cg = Path('/sys/fs/cgroup') / group.lstrip('/')
    limits = {k: (cg / k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert limits == {'cpu.max': '200000 100000', 'memory.max': str(24 * 2**30), 'memory.swap.max': '0'}
    assert all(os.environ.get(k) == '1' for k in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'])
    jobs = json.loads(Path(plan['jobs']).read_text()); assert len(jobs) == 1620
    root = Path(plan['output']); root.mkdir(exist_ok=True)
    lock = (root / 'stage.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (root / 'receipt.json').exists(), 'Completed startup stage cannot restart'
    marker = root / 'stage_plan.json'
    if marker.exists(): assert json.loads(marker.read_text()) == plan
    else: write_json(marker, plan)
    (root / 'chains').mkdir(exist_ok=True); results = []; started = time.monotonic()
    def run(job):
        chain = job['chain']; cp = root / 'chains' / (chain['chain_id'] + '.json')
        assert shutil.disk_usage(root).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
        if cp.exists():
            old = json.loads(cp.read_text()); assert old['plan_sha256'] == digest
            reconstructed = inspect(job, Path(old['attempt_receipt']), digest)
            assert old == reconstructed, 'Startup checkpoint changed'; return old
        folder = root / 'attempts' / chain['chain_id']
        receipts = sorted(folder.glob('attempt-*/receipt.json'))
        if receipts:
            assert len(receipts) == 1, 'Never silently choose between startup attempts'
            receipt = receipts[0]
        else: receipt = run_attempt(folder, job['config'])
        row = inspect(job, receipt, digest); write_json(cp, row); return row
    with ThreadPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(run, job) for job in jobs]):
            row = future.result(); results.append(row)
            print('reference_native_startups', len(results), '/1620', row['chain_id'], row['status'], flush=True)
    verify(plan); assert sha(a.plan) == digest
    results.sort(key=lambda x: x['chain_id']); write_json(root / 'dispositions.json', results)
    bindings = dict(plan['pins']); bindings[str(a.plan)] = digest
    artifacts = {str(p.relative_to(root)): sha(p) for p in [marker, root / 'dispositions.json',
                 *sorted((root / 'chains').glob('*.json'))]}
    receipt = dict(status='complete_full_reference_startup_dispositions_pending_readback',
        plan_sha256=digest, **summarize(results), source_hashes=bindings, artifacts=artifacts,
        elapsed_seconds=time.monotonic() - started, scientific_eligibility=False, scope=plan['scope'])
    write_json(root / 'receipt.json', receipt)
    print(json.dumps({k: receipt[k] for k in SUMMARY_FIELDS}), flush=True)


if __name__ == '__main__': main()
