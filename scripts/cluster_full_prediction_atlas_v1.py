#!/usr/bin/env python3
"""Cluster a fully validated prediction atlas and verify complete membership.

This runner is deliberately unusable until a launch-specific plan binds completed
database, coordinate-profile, and PAE readback receipts.  It creates a candidate
structural partition only; it makes no homology or evolutionary claim.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import time


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def require_receipt(spec):
    path = Path(spec['path'])
    if sha(path) != spec['sha256']:
        raise ValueError('Receipt hash differs: ' + str(path))
    receipt = json.loads(path.read_text())
    if receipt.get('status') != spec['status']:
        raise ValueError('Unexpected receipt status: ' + str(path))
    return receipt


def partition_counts(path, expected):
    assigned, representatives, counts = set(), set(), Counter()
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 2:
                raise ValueError(f'Invalid cluster row {number}')
            representative, member = fields
            if representative not in expected or member not in expected:
                raise ValueError(f'Unknown cluster member at row {number}')
            if member in assigned:
                raise ValueError(f'Repeated cluster member at row {number}')
            assigned.add(member)
            representatives.add(representative)
            counts[representative] += 1
    if assigned != expected:
        raise ValueError('Cluster partition does not cover every database model')
    if not representatives.issubset(expected):
        raise ValueError('Unknown cluster representative')
    return counts


def self_test():
    expected = {'a', 'b', 'c'}
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'members.tsv'
        path.write_text('a\ta\na\tb\nc\tc\n')
        counts = partition_counts(path, expected)
        assert counts == {'a': 2, 'c': 1}
        path.write_text('a\ta\na\tb\na\tb\n')
        try:
            partition_counts(path, expected)
        except ValueError:
            pass
        else:
            raise AssertionError('Duplicate member accepted')
    return ['complete_partition_membership', 'duplicate_member_rejected']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    checks = self_test()
    if args.self_test:
        print(json.dumps({'status': 'passed_full_atlas_clustering_runner_controls', 'checks': checks}, indent=2))
        return
    if args.plan is None:
        parser.error('--plan is required unless --self-test is used')

    plan = json.loads(args.plan.read_text())
    for path, digest in plan['pins'].items():
        if sha(path) != digest:
            raise ValueError('Pinned input changed: ' + path)
    receipts = {name: require_receipt(spec) for name, spec in plan['required_receipts'].items()}
    database = Path(plan['database'])
    lookup = database / plan['lookup_name']
    expected = set()
    with lookup.open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 3:
                raise ValueError(f'Invalid database lookup row {number}')
            name = fields[1]
            if name in expected:
                raise ValueError('Repeated database source alias')
            expected.add(name)
    if len(expected) != plan['expected_models']:
        raise ValueError('Database model count differs from frozen plan')
    for name, digest in plan['database_artifacts'].items():
        if sha(database / name) != digest:
            raise ValueError('Database artifact changed: ' + name)

    output = Path(plan['output'])
    if output.exists():
        raise FileExistsError('Fresh output required: ' + str(output))
    resource = plan['resources']
    if shutil.disk_usage(output.parent).free < resource['minimum_free_disk_gib'] * 2**30:
        raise RuntimeError('Insufficient free disk for full clustering')
    output.mkdir(parents=True)
    state = output / 'state.json'
    prefix, cluster = str(database / plan['database_prefix']), str(output / 'clusters')
    commands = [
        [plan['foldseek'], 'cluster', prefix, cluster, str(output / 'tmp'), *plan['cluster_arguments']],
        [plan['foldseek'], 'createtsv', prefix, prefix, cluster, str(output / 'cluster_members.tsv'),
         '--threads', str(resource['cpu'])],
    ]
    try:
        for stage, command in enumerate(commands):
            state.write_text(json.dumps({'status': 'running_native_full_atlas_clustering',
                                         'stage': stage, 'command': command}, indent=2) + '\n')
            with (output / f'stage_{stage}.log').open('w') as log:
                child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                         env=dict(os.environ, CUDA_VISIBLE_DEVICES=''), start_new_session=True)
                while child.poll() is None:
                    if shutil.disk_usage(output).free < resource['emergency_free_disk_gib'] * 2**30:
                        os.killpg(child.pid, signal.SIGTERM)
                        raise RuntimeError('Emergency disk reserve reached')
                    time.sleep(10)
                if child.returncode:
                    raise RuntimeError(f'Native Foldseek stage {stage} exited {child.returncode}')
        counts = partition_counts(output / 'cluster_members.tsv', expected)
        with (output / 'cluster_sizes.tsv').open('w') as handle:
            handle.write('representative\tmodels\n')
            for representative, count in sorted(counts.items()):
                handle.write(f'{representative}\t{count}\n')
        for name, digest in plan['database_artifacts'].items():
            if sha(database / name) != digest:
                raise ValueError('Database changed during native clustering: ' + name)
        result = {
            'status': 'complete_full_atlas_candidate_structural_partition_pending_independent_readback',
            'models': len(expected), 'clusters': len(counts),
            'singleton_clusters': sum(count == 1 for count in counts.values()),
            'largest_cluster_models': max(counts.values()), 'plan_sha256': sha(args.plan),
            'required_receipt_sha256': {name: sha(spec['path']) for name, spec in plan['required_receipts'].items()},
            'artifacts': {p.name: sha(p) for p in output.iterdir() if p.is_file() and p.name != 'state.json'},
            'scientific_eligibility': False,
            'scope': 'Complete native candidate structural partition with exact database-membership checks. '
                     'It does not establish structural homology, orthology, confidence calibration, '
                     'evolutionary events, or a biological result.',
        }
        (output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
        state.write_text(json.dumps({'status': 'complete', 'models': len(expected), 'clusters': len(counts)}, indent=2) + '\n')
        print(json.dumps(result, indent=2))
    except BaseException:
        state.write_text(json.dumps({'status': 'failed_or_interrupted_without_accepted_receipt'}, indent=2) + '\n')
        raise


if __name__ == '__main__':
    main()
