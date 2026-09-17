#!/usr/bin/env python3
"""Exhaustively read back all empirical overlap pairs without using the sweep kernel."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def exhaustive(hits, nested):
    result = {}
    order = sorted(hits, key=lambda h: (h['envelope_start'], h['envelope_end'], h['hit_id']))
    for a, b in itertools.combinations(order, 2):
        left = max(a['envelope_start'], b['envelope_start'])
        right = min(a['envelope_end'], b['envelope_end'])
        if right < left:
            continue
        al = max(a['alignment_start'], b['alignment_start'])
        ar = min(a['alignment_end'], b['alignment_end'])
        ac = a['alignment_start'] <= b['alignment_start'] <= b['alignment_end'] <= a['alignment_end']
        bc = b['alignment_start'] <= a['alignment_start'] <= a['alignment_end'] <= b['alignment_end']
        an = (a['pfam_accession'], b['pfam_accession']) in nested
        bn = (b['pfam_accession'], a['pfam_accession']) in nested
        result[(a['hit_id'], b['hit_id'])] = dict(alignment_overlap=max(0, ar-al+1),
            envelope_overlap=right-left+1, same_family=int(a['pfam_accession']==b['pfam_accession']),
            same_clan=int(bool(a['pfam_clan']) and a['pfam_clan']==b['pfam_clan']),
            curated_a_contains_b=int(an), curated_b_contains_a=int(bn),
            alignment_a_contains_b=int(ac), alignment_b_contains_a=int(bc),
            curated_and_coordinate_nesting=int((an and ac) or (bn and bc)))
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    plan = json.loads(a.plan.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    root = Path(plan['output'])
    receipt = json.loads((root/'receipt.json').read_text())
    assert receipt['plan_sha256'] == sha(a.plan)
    assert receipt['status'] == 'complete_full_overlap_inventory_requires_independent_readback'
    for name, digest in receipt['artifacts'].items():
        assert sha(root/name) == digest
    database = Path(plan['database_root'])
    dr = json.loads((database/'receipt.json').read_text())
    assert sha(database/'domains.sqlite') == dr['artifacts']['domains.sqlite']
    with Path(plan['nested_table']).open() as handle:
        nested = {(r['containing_accession'],r['nested_accession']) for r in csv.DictReader(handle,delimiter='\t')}
    conn = sqlite3.connect('file:'+str((database/'domains.sqlite').resolve())+'?mode=ro',uri=True)
    conn.row_factory = sqlite3.Row
    fields = ['hit_id','pfam_accession','pfam_clan','pfam_type','alignment_start','alignment_end','envelope_start','envelope_end','hmm_coverage']
    sql = 'SELECT q.sequence_id,q.search_partition,'+','.join('h.'+k for k in fields)+' FROM queries q LEFT JOIN hits h ON h.sequence_id=q.sequence_id ORDER BY q.sequence_id'
    totals = Counter()
    with gzip.open(root/'query_annotation_inventory.tsv.gz','rt') as qh, gzip.open(root/'overlap_pairs.tsv.gz','rt') as ph:
        inventory = iter(csv.DictReader(qh,delimiter='\t'))
        pair_rows = iter(csv.DictReader(ph,delimiter='\t'))
        current = next(pair_rows, None)
        for sequence, group in itertools.groupby(conn.execute(sql),key=lambda r:r['sequence_id']):
            rows = list(group)
            partition = rows[0]['search_partition']
            hits = [dict(r) for r in rows if r['hit_id'] is not None]
            for h in hits:
                for k in ['alignment_start','alignment_end','envelope_start','envelope_end']:
                    h[k] = int(h[k])
            expected = exhaustive(hits,nested)
            actual = {}
            while current is not None and current['sequence_id']==sequence:
                assert current['search_partition']==partition
                key = (current['hit_a'],current['hit_b'])
                assert key not in actual
                actual[key] = {k:int(v) for k,v in current.items() if k not in ['sequence_id','search_partition','hit_a','hit_b']}
                current = next(pair_rows,None)
            assert actual==expected, sequence
            if current is not None:
                assert current['sequence_id']>sequence
            item = next(inventory,None)
            assert item is not None and item['sequence_id']==sequence and item['search_partition']==partition
            align_count = sum(v['alignment_overlap']>0 for v in expected.values())
            assert int(item['hits'])==len(hits)
            assert int(item['alignment_overlap_pairs'])==align_count
            assert int(item['envelope_overlap_pairs'])==len(expected)
            assert int(item['partial_hmm_hits_below_070'])==sum(float(h['hmm_coverage'])<.7 for h in hits)
            ordered = sorted(hits,key=lambda h:(h['alignment_start'],h['alignment_end'],h['hit_id']))
            assert json.loads(item['ordered_raw_annotations_json'])==[[h['hit_id'],h['pfam_accession'],h['pfam_type'],h['alignment_start'],h['alignment_end']] for h in ordered]
            status = 'no_GA_hit_not_proven_absence' if not hits else 'overlaps_require_review' if expected else 'nonoverlapping_raw_annotations_not_validated_architecture'
            assert item['status']==status
            totals.update(queries=1,hits=len(hits),alignment_overlap_pairs=align_count,envelope_overlap_pairs=len(expected))
            totals[status]+=1
            totals['same_clan_pairs']+=sum(v['same_clan'] for v in expected.values())
            totals['curated_and_coordinate_nesting_pairs']+=sum(v['curated_and_coordinate_nesting'] for v in expected.values())
            if totals['queries']%100000==0:
                print(totals['queries'],'queries independently checked',flush=True)
        assert next(inventory,None) is None and current is None
    conn.close()
    assert totals==Counter(receipt['totals'])
    assert totals['queries']==dr['queries'] and totals['hits']==dr['total_hit_rows']
    result = dict(status='passed_complete_exhaustive_overlap_readback',totals=dict(totals),
                  source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(Path(__file__)),
                  scope='Every query inventory and every envelope-overlap pair checked against exhaustive within-query combinations, including all coordinate, clan and directed nesting flags. Does not resolve competing annotations or establish biological architectures.')
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
