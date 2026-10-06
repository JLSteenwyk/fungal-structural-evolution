#!/usr/bin/env python3
"""Independently replay a full-atlas candidate structural cluster partition."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import tempfile


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def lookup_members(path, expected_models):
    members = set()
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 3 or not fields[1]:
                raise ValueError(f'Malformed lookup record {number}')
            if fields[1] in members:
                raise ValueError(f'Repeated lookup identifier at row {number}')
            members.add(fields[1])
    if len(members) != expected_models:
        raise ValueError('Lookup count differs from expected full atlas')
    return members


def replay_partition(path, expected):
    assigned, self_members, sizes = set(), set(), Counter()
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 2 or not all(fields):
                raise ValueError(f'Malformed partition record {number}')
            representative, member = fields
            if representative not in expected or member not in expected:
                raise ValueError(f'Unknown identifier at partition record {number}')
            if member in assigned:
                raise ValueError(f'Repeated member at partition record {number}')
            assigned.add(member)
            sizes[representative] += 1
            if representative == member:
                self_members.add(member)
    if assigned != expected:
        raise ValueError('Partition membership differs from database lookup')
    if set(sizes) != self_members:
        raise ValueError('Representative self-membership differs from cluster set')
    return sizes


def controls():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory); lookup = root / 'lookup'; partition = root / 'partition'
        lookup.write_text('1\ta\t0\n2\tb\t1\n3\tc\t2\n')
        expected = lookup_members(lookup, 3)
        partition.write_text('a\ta\na\tb\nc\tc\n')
        assert replay_partition(partition, expected) == {'a': 2, 'c': 1}
        partition.write_text('a\ta\na\tb\nc\tb\n')
        try:
            replay_partition(partition, expected)
        except ValueError:
            pass
        else:
            raise AssertionError('Repeated member accepted')
        partition.write_text('a\ta\na\tb\nc\tc\n')
    return ['independent_lookup_parse', 'complete_partition_replay', 'duplicate_member_rejected']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    checks = controls()
    if args.self_test:
        print(json.dumps({'status': 'passed_full_atlas_cluster_readback_controls', 'checks': checks}, indent=2))
        return
    if args.plan is None:
        parser.error('--plan is required unless --self-test is used')
    plan = json.loads(args.plan.read_text())
    for path, expected in plan['pins'].items():
        if digest(path) != expected:
            raise ValueError('Pinned input changed: ' + path)
    producer_receipt = Path(plan['producer_receipt'])
    producer = json.loads(producer_receipt.read_text())
    if producer.get('status') != 'complete_full_atlas_candidate_structural_partition_pending_independent_readback':
        raise ValueError('Completed full-atlas clustering producer receipt required')
    if digest(producer_receipt) != plan['producer_receipt_sha256']:
        raise ValueError('Producer receipt hash differs')
    output = Path(plan['producer_output'])
    for name, expected in producer['artifacts'].items():
        if digest(output / name) != expected:
            raise ValueError('Producer artifact changed: ' + name)
    database = Path(plan['database'])
    for name, expected in plan['database_artifacts'].items():
        if digest(database / name) != expected:
            raise ValueError('Database artifact changed: ' + name)
    expected = lookup_members(database / plan['lookup_name'], plan['expected_models'])
    sizes = replay_partition(output / 'cluster_members.tsv', expected)
    if producer['models'] != len(expected) or producer['clusters'] != len(sizes):
        raise ValueError('Producer summary differs from independent replay')
    if producer['singleton_clusters'] != sum(n == 1 for n in sizes.values()):
        raise ValueError('Producer singleton count differs from independent replay')
    if producer['largest_cluster_models'] != max(sizes.values()):
        raise ValueError('Producer largest-cluster count differs from independent replay')
    for name, expected_hash in plan['database_artifacts'].items():
        if digest(database / name) != expected_hash:
            raise ValueError('Database changed during reader replay: ' + name)
    result = dict(status='passed_independent_full_atlas_candidate_structural_partition_readback',
                  models=len(expected), clusters=len(sizes),
                  singleton_clusters=sum(n == 1 for n in sizes.values()),
                  largest_cluster_models=max(sizes.values()), checks=checks,
                  producer_receipt_sha256=digest(producer_receipt), plan_sha256=digest(args.plan),
                  scientific_eligibility=False,
                  scope='Independent lookup and complete partition-membership replay. It does not validate '\
                        'alignment thresholds, confidence/PAE calibration, structural homology, orthology, '\
                        'or evolutionary events.')
    receipt = Path(plan['output'])
    if receipt.exists():
        raise FileExistsError('Fresh reader receipt required')
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
