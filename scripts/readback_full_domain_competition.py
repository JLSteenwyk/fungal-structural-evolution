#!/usr/bin/env python3
"""Audit all four competition policies with independently constructed pair graphs."""
import argparse
from collections import Counter
import csv
from decimal import Decimal
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def independent_policy(hits,nested,coord,rank):
    """Construct all pair relationships before solving ordered competition."""
    n=len(hits)
    ids=[h['hit_id'] for h in hits]
    metrics=[(Decimal(h['independent_evalue']),-Decimal(h['domain_score'])) for h in hits]
    scores=metrics if rank=='evalue' else [(b,e) for e,b in metrics]
    priority=sorted(range(n),key=lambda i:(*scores[i],ids[i]))
    graph=[set() for _ in hits]
    pair_info=[]
    ties=[]
    for i,j in itertools.combinations(range(n),2):
        a,b=hits[i],hits[j]
        intersects=not (int(a[coord+'_end'])<int(b[coord+'_start']) or int(b[coord+'_end'])<int(a[coord+'_start']))
        if not intersects:
            continue
        ab=((a['pfam_accession'],b['pfam_accession']) in nested
            and int(a['alignment_start'])<int(b['alignment_start'])
            and int(b['alignment_end'])<int(a['alignment_end']))
        ba=((b['pfam_accession'],a['pfam_accession']) in nested
            and int(b['alignment_start'])<int(a['alignment_start'])
            and int(a['alignment_end'])<int(b['alignment_end']))
        is_nested=ab or ba
        same=bool(a['pfam_clan']) and a['pfam_clan']==b['pfam_clan']
        pair_info.append((i,j,is_nested))
        if same and not is_nested:
            graph[i].add(j)
            graph[j].add(i)
            if scores[i][0]==scores[j][0]:
                ties.append(sorted([ids[i],ids[j]]))
    kept=set()
    decisions=[]
    for i in priority:
        blockers=[j for j in priority if j in kept and j in graph[i]]
        if not blockers:
            kept.add(i)
        decisions.append(dict(hit_id=ids[i],disposition='suppressed_clan_competitor' if blockers else 'retained_candidate',
                              competing_retained_hits=[ids[j] for j in blockers]))
    ordered=sorted(kept,key=lambda i:(int(hits[i]['alignment_start']),int(hits[i]['alignment_end']),ids[i]))
    unresolved=[]
    nesting=[]
    for i,j,nested_pair in pair_info:
        if i in kept and j in kept:
            (nesting if nested_pair else unresolved).append(sorted([ids[i],ids[j]]))
    return dict(retained_hits=[ids[i] for i in ordered],decisions=sorted(decisions,key=lambda d:d['hit_id']),
                primary_rank_ties=sorted(ties),unresolved_overlap_pairs=sorted(unresolved),
                candidate_nested_pairs=sorted(nesting),coordinates=coord,ranking=rank)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    plan=json.loads(a.plan.read_text())
    for path,digest in plan['pins'].items():
        assert sha(path)==digest,path
    root=Path(plan['output'])
    receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='completed_full_competition_requires_independent_readback'
    assert receipt['plan_sha256']==sha(a.plan)
    for name,digest in receipt['artifacts'].items():
        assert sha(root/name)==digest
    source=Path(plan['database_root'])
    dbreceipt=json.loads((source/'receipt.json').read_text())
    assert sha(source/'domains.sqlite')==dbreceipt['artifacts']['domains.sqlite']
    with Path(plan['nested_table']).open() as handle:
        nested={(r['containing_accession'],r['nested_accession']) for r in csv.DictReader(handle,delimiter='\t')}
    conn=sqlite3.connect('file:'+str((source/'domains.sqlite').resolve())+'?mode=ro',uri=True)
    conn.row_factory=sqlite3.Row
    sql='SELECT q.sequence_id AS query_id,q.search_partition AS query_partition,h.* FROM queries q LEFT JOIN hits h ON h.sequence_id=q.sequence_id ORDER BY q.sequence_id'
    totals=Counter()
    policy_totals={k:Counter() for k in receipt['policies']}
    assert set(policy_totals)=={c+'_'+r for c in ['alignment','envelope'] for r in ['evalue','bitscore']}
    csv.field_size_limit(256*1024*1024)
    with gzip.open(root/'query_competition.tsv.gz','rt') as handle:
        output=csv.DictReader(handle,delimiter='\t')
        for query,group in itertools.groupby(conn.execute(sql),key=lambda r:r['query_id']):
            rows=list(group)
            hits=[dict(r) for r in rows if r['hit_id'] is not None]
            observed=next(output)
            assert observed['sequence_id']==query
            assert observed['search_partition']==rows[0]['query_partition']
            assert int(observed['raw_hits'])==len(hits)
            partial=sum(Decimal(h['hmm_coverage'])<Decimal('0.70') for h in hits)
            assert int(observed['raw_partial_hmm_hits_below_070'])==partial
            assert observed['status']==('candidate_annotations_require_validation' if hits else 'no_GA_hit_not_proven_absence')
            actual=json.loads(observed['policies_json'])
            assert set(actual)==set(policy_totals)
            sets=[]
            for key in policy_totals:
                coord,rank=key.split('_')
                expected=independent_policy(hits,nested,coord,rank)
                assert actual[key]==expected,(query,key)
                pc=policy_totals[key]
                pc['retained_hits']+=len(expected['retained_hits'])
                pc['suppressed_hits']+=len(hits)-len(expected['retained_hits'])
                pc['queries_with_primary_rank_ties']+=bool(expected['primary_rank_ties'])
                pc['queries_with_unresolved_overlaps']+=bool(expected['unresolved_overlap_pairs'])
                pc['queries_with_candidate_nesting']+=bool(expected['candidate_nested_pairs'])
                sets.append(frozenset(expected['retained_hits']))
            disagree=any(s!=sets[0] for s in sets)
            assert int(observed['retained_sets_agree'])==int(not disagree)
            totals.update(queries=1,hits=len(hits),partial_hmm_hits_below_070=partial,
                          queries_with_policy_disagreement=int(disagree),queries_without_hits=int(not hits))
            if totals['queries']%100000==0:
                print(totals['queries'],'queries independently checked',flush=True)
        assert next(output,None) is None
    conn.close()
    assert dict(totals)==receipt['totals']
    assert {k:dict(v) for k,v in policy_totals.items()}==receipt['policies']
    assert totals['queries']==dbreceipt['queries'] and totals['hits']==dbreceipt['total_hit_rows']
    result=dict(status='passed_complete_independent_domain_competition_readback',totals=dict(totals),
                policies={k:dict(v) for k,v in policy_totals.items()},plan_sha256=sha(a.plan),
                production_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(Path(__file__)),
                scope='Every query and every four-policy disposition, retained ordering, blocker list, primary tie, residual overlap and nested candidate compared to an independently constructed exhaustive pair graph. All summaries matched. Not biological architecture or gain/loss validation.')
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    main()
