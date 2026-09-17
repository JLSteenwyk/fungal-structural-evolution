#!/usr/bin/env python3
"""Map all hit overlaps, retaining clan and curated nesting annotations without filtering."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import sqlite3


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def overlap(a, b, start, end):
    return max(0, min(a[end], b[end]) - max(a[start], b[start]) + 1)


def describe_pair(a, b, nested):
    aln = overlap(a, b, 'alignment_start', 'alignment_end')
    env = overlap(a, b, 'envelope_start', 'envelope_end')
    forward = (a['pfam_accession'], b['pfam_accession']) in nested
    reverse = (b['pfam_accession'], a['pfam_accession']) in nested
    contains_ab = a['alignment_start'] <= b['alignment_start'] and a['alignment_end'] >= b['alignment_end']
    contains_ba = b['alignment_start'] <= a['alignment_start'] and b['alignment_end'] >= a['alignment_end']
    return dict(hit_a=a['hit_id'], hit_b=b['hit_id'], alignment_overlap=aln,
                envelope_overlap=env, same_family=int(a['pfam_accession'] == b['pfam_accession']),
                same_clan=int(bool(a['pfam_clan']) and a['pfam_clan'] == b['pfam_clan']),
                curated_a_contains_b=int(forward), curated_b_contains_a=int(reverse),
                alignment_a_contains_b=int(contains_ab), alignment_b_contains_a=int(contains_ba),
                curated_and_coordinate_nesting=int((forward and contains_ab) or (reverse and contains_ba)))


def pairs(hits, nested):
    active = []
    for b in sorted(hits, key=lambda h: (h['envelope_start'], h['envelope_end'], h['hit_id'])):
        active = [a for a in active if a['envelope_end'] >= b['envelope_start']]
        for a in active:
            yield describe_pair(a, b, nested)
        active.append(b)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    a = ap.parse_args()
    plan = json.loads(a.plan.read_text())
    for path, digest in plan['pins'].items():
        if sha(path) != digest:
            raise ValueError('Changed pinned file: ' + path)
    source = Path(plan['database_root'])
    receipt = json.loads((source / 'receipt.json').read_text())
    audit = json.loads((source / 'readback.json').read_text())
    assert audit['status'] == 'passed_complete_domain_database_source_readback'
    assert audit['database_receipt_sha256'] == sha(source / 'receipt.json')
    dbpath = source / 'domains.sqlite'
    assert sha(dbpath) == receipt['artifacts']['domains.sqlite']
    with Path(plan['nested_table']).open() as handle:
        nested = {(r['containing_accession'], r['nested_accession']) for r in csv.DictReader(handle, delimiter='\t')}
    output = Path(plan['output'])
    if output.exists():
        raise FileExistsError(output)
    assert shutil.disk_usage(output.parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    output.mkdir(parents=True)
    conn = sqlite3.connect('file:' + str(dbpath.resolve()) + '?mode=ro', uri=True)
    conn.row_factory = sqlite3.Row
    fields = ['hit_id','pfam_accession','pfam_clan','pfam_type','alignment_start','alignment_end','envelope_start','envelope_end','hmm_coverage']
    sql = 'SELECT q.sequence_id,q.search_partition,' + ','.join('h.'+f for f in fields) + ' FROM queries q LEFT JOIN hits h ON h.sequence_id=q.sequence_id ORDER BY q.sequence_id'
    cursor = conn.execute(sql)
    pair_path, query_path = output / 'overlap_pairs.tsv.gz', output / 'query_annotation_inventory.tsv.gz'
    totals = Counter()
    pair_fields = ['sequence_id','search_partition','hit_a','hit_b','alignment_overlap','envelope_overlap','same_family','same_clan','curated_a_contains_b','curated_b_contains_a','alignment_a_contains_b','alignment_b_contains_a','curated_and_coordinate_nesting']
    query_fields = ['sequence_id','search_partition','hits','alignment_overlap_pairs','envelope_overlap_pairs','partial_hmm_hits_below_070','status','ordered_raw_annotations_json']
    with gzip.open(pair_path, 'wt') as ph, gzip.open(query_path, 'wt') as qh:
        pw = csv.DictWriter(ph, fieldnames=pair_fields, delimiter='\t', lineterminator='\n')
        qw = csv.DictWriter(qh, fieldnames=query_fields, delimiter='\t', lineterminator='\n')
        pw.writeheader()
        qw.writeheader()
        for sequence, group in itertools.groupby(cursor, key=lambda r: r['sequence_id']):
            rows = list(group)
            partition = rows[0]['search_partition']
            hits = [dict(r) for r in rows if r['hit_id'] is not None]
            for h in hits:
                for field in ['alignment_start','alignment_end','envelope_start','envelope_end']:
                    h[field] = int(h[field])
                assert 1 <= h['envelope_start'] <= h['alignment_start'] <= h['alignment_end'] <= h['envelope_end']
            aln_pairs = env_pairs = 0
            for pair in pairs(hits, nested):
                env_pairs += 1
                aln_pairs += pair['alignment_overlap'] > 0
                totals['same_clan_pairs'] += pair['same_clan']
                totals['curated_and_coordinate_nesting_pairs'] += pair['curated_and_coordinate_nesting']
                pw.writerow(dict(sequence_id=sequence, search_partition=partition, **pair))
            status = 'no_GA_hit_not_proven_absence' if not hits else 'overlaps_require_review' if env_pairs else 'nonoverlapping_raw_annotations_not_validated_architecture'
            ordered = sorted(hits, key=lambda h: (h['alignment_start'], h['alignment_end'], h['hit_id']))
            annotations = [[h['hit_id'],h['pfam_accession'],h['pfam_type'],h['alignment_start'],h['alignment_end']] for h in ordered]
            qw.writerow(dict(sequence_id=sequence, search_partition=partition, hits=len(hits),
                             alignment_overlap_pairs=aln_pairs, envelope_overlap_pairs=env_pairs,
                             partial_hmm_hits_below_070=sum(float(h['hmm_coverage']) < .7 for h in hits),
                             status=status, ordered_raw_annotations_json=json.dumps(annotations, separators=(',', ':'))))
            totals.update(queries=1, hits=len(hits), alignment_overlap_pairs=aln_pairs, envelope_overlap_pairs=env_pairs)
            totals[status] += 1
            if totals['queries'] % 100000 == 0:
                print(totals['queries'], 'queries', totals['hits'], 'hits', flush=True)
                if sum(p.stat().st_size for p in output.iterdir()) > plan['resources']['output_allowance_gib'] * 2**30:
                    raise ValueError('Output allowance exceeded')
    conn.close()
    assert totals['queries'] == receipt['queries'] and totals['hits'] == receipt['total_hit_rows']
    result = dict(status='complete_full_overlap_inventory_requires_independent_readback', totals=dict(totals),
                  plan_sha256=sha(a.plan), script_sha256=sha(Path(__file__)),
                  artifacts={p.name: sha(p) for p in [pair_path, query_path]},
                  scope='Every query and source hit retained; all inclusive envelope-overlap pairs enumerated with alignment overlap, clan and curated nesting flags. No competing hit removed, no automatic nested-domain acceptance, no architecture or gain/loss inference.')
    (output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
