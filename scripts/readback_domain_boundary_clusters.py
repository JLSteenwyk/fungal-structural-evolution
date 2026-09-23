#!/usr/bin/env python3
"""Independently reconstruct every domain boundary/cluster output row."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import time
import psutil


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        yield from csv.DictReader(f, delimiter='\t')


def audit_tables(manifest, lookup, members, output):
    intervals = {}
    for r in read_rows(manifest / 'intervals.tsv'):
        key = r['interval_id']
        if key in intervals:
            raise ValueError('Repeated source interval')
        intervals[key] = r['model_key']
    expected = set()
    with lookup.open() as f:
        for line in f:
            _, key, _ = line.rstrip('\n').split('\t')
            if key in expected or key not in intervals:
                raise ValueError('Invalid lookup interval')
            expected.add(key)
    assignment = {}
    with members.open() as f:
        for line in f:
            rep, member = line.rstrip('\n').split('\t')
            if member in assignment or member not in expected or rep not in expected:
                raise ValueError('Invalid partition member')
            assignment[member] = rep
    if set(assignment) != expected or any(assignment[r] != r for r in assignment.values()):
        raise ValueError('Incomplete partition or non-self representative')
    pairs = {}; linked = set(); links = 0
    for r in read_rows(manifest / 'boundary_links.tsv.gz'):
        key = r['model_key'], r['hit_id']; interval = r['interval_id']
        if interval not in intervals or intervals[interval] != key[0]:
            raise ValueError('Boundary refers to another model or unknown interval')
        slot = {'alignment': 0, 'envelope': 1}.get(r['boundary'])
        bounds = pairs.setdefault(key, [None, None])
        if slot is None or bounds[slot] is not None:
            raise ValueError('Unknown or duplicate boundary')
        bounds[slot] = interval; linked.add(interval); links += 1
    if linked != set(intervals) or any(None in bounds for bounds in pairs.values()):
        raise ValueError('Incomplete boundary grid')
    seen = set(); counts = Counter()
    for r in read_rows(output):
        key = r['model_key'], r['hit_id']
        if key in seen or key not in pairs:
            raise ValueError('Repeated or unexpected result pair')
        a, e = pairs[key]; ac = assignment.get(a, ''); ec = assignment.get(e, '')
        # Independent decision order: classify by interval presence before identities.
        presence = int(a in assignment) + 2 * int(e in assignment)
        if presence < 3:
            category = {0: 'both_unclustered', 1: 'envelope_unclustered',
                        2: 'alignment_unclustered'}[presence]
        elif a == e:
            category = 'identical_interval'
        elif ac != ec:
            category = 'distinct_intervals_different_clusters'
        else:
            category = 'distinct_intervals_same_cluster'
        reconstructed = dict(model_key=key[0], hit_id=key[1], alignment_interval=a,
                             envelope_interval=e, alignment_cluster=ac,
                             envelope_cluster=ec, disposition=category)
        if r != reconstructed:
            raise ValueError('Output field differs from source reconstruction: '+repr(key))
        seen.add(key); counts[category] += 1
    if seen != set(pairs):
        raise ValueError('Missing output pairs')
    return dict(candidate_pairs=len(pairs), intervals=len(intervals),
                clustered_intervals=len(expected), unclustered_intervals=len(intervals)-len(expected),
                boundary_links=links, dispositions=dict(counts))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', required=True, type=Path)
    args = ap.parse_args(); plan = json.loads(args.plan.read_text()); ph = sha(args.plan)
    pins = dict(plan['pins']); pins[str(args.plan)] = ph
    def verify():
        for path, digest in pins.items():
            if sha(path) != digest:
                raise ValueError('Changed pinned source: '+path)
    verify(); out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    state = out / 'state.json'
    state.write_text(json.dumps({'status': 'waiting_for_boundary_comparison'})+'\n')
    dependency = plan['predecessor']
    while True:
        try:
            proc = psutil.Process(dependency['pid'])
            live = proc.create_time() == dependency['create_time'] and proc.status() != psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:
            live = False
        if not live:
            break
        time.sleep(20)
    verify()
    if shutil.disk_usage(out).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
        raise ValueError('Insufficient free disk')
    producer = json.loads(Path(plan['producer_plan']).read_text())
    roots = {k: Path(producer[k]) for k in ['output', 'clusters', 'database', 'manifest']}
    receipts = {}
    for key, root in roots.items():
        path = root / 'receipt.json'; pins[str(path)] = sha(path)
        receipts[key] = json.loads(path.read_text())
    r, cr, dr, mr = (receipts[k] for k in ['output', 'clusters', 'database', 'manifest'])
    if (r['status'] != 'complete_all_domain_boundary_cluster_dispositions'
            or r['plan_sha256'] != sha(plan['producer_plan'])
            or r['cluster_receipt_sha256'] != sha(roots['clusters'] / 'receipt.json')
            or r['database_receipt_sha256'] != sha(roots['database'] / 'receipt.json')
            or cr['status'] != 'complete_domain_candidate_partition_membership_readback'
            or cr['plan_sha256'] != sha(producer['cluster_plan'])
            or cr['database_receipt_sha256'] != sha(roots['database'] / 'receipt.json')
            or dr['status'] != 'complete_domain_foldseek_database_with_full_sequence_and_coordinate_readback'
            or dr['plan_sha256'] != sha(producer['database_plan'])
            or mr['status'] != 'complete_all_candidate_domain_extraction_manifest'):
        raise ValueError('Incomplete or incorrectly bound source receipt')
    lookup = roots['database'] / 'domains.lookup'; members = roots['clusters'] / 'cluster_members.tsv'
    table = roots['output'] / 'boundary_cluster_dispositions.tsv.gz'
    for path, digest in [(lookup, dr['artifacts'][lookup.name]),
                         (members, cr['artifacts'][members.name]), (table, r['artifacts'][table.name])]:
        pins[str(path)] = digest
    verify(); state.write_text(json.dumps({'status': 'reconstructing_all_rows'})+'\n')
    counts = audit_tables(roots['manifest'], lookup, members, table)
    if (any(r[k] != v for k, v in counts.items() if k != 'boundary_links')
            or counts['candidate_pairs'] != mr['candidate_model_hit_pairs']
            or counts['intervals'] != mr['unique_intervals']
            or counts['boundary_links'] != mr['boundary_links']
            or counts['clustered_intervals'] != cr['intervals']
            or counts['clustered_intervals'] != dr['models']
            or counts['unclustered_intervals'] != cr['excluded_rejected_intervals']):
        raise ValueError('Full scope or summary counts differ')
    verify()
    result = dict(status='passed_full_domain_boundary_cluster_output_readback', **counts,
                  plan_sha256=ph, source_hashes=pins,
                  scope='Every result field reconstructed from original interval identities, paired boundaries, database lookup and complete partition without producer comparison functions. Checks within-partition assignment sensitivity only; no independent clustering, alignment-threshold, homology, biological boundary or evolutionary validation.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    state.write_text(json.dumps({'status': result['status']})+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
