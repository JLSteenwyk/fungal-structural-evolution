#!/usr/bin/env python3
"""Stream all domain queries through four explicit clan-competition policies."""
import argparse
from collections import Counter
import csv
from decimal import Decimal
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import sqlite3
import time
from domain_clan_competition import compete


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan',type=Path,required=True)
    a=ap.parse_args()
    plan=json.loads(a.plan.read_text())
    for path,digest in plan['pins'].items():
        if sha(path)!=digest:
            raise ValueError('Changed pin: '+path)
    root=Path(plan['database_root'])
    receipt=json.loads((root/'receipt.json').read_text())
    audit=json.loads((root/'readback.json').read_text())
    if (audit['status']!='passed_complete_domain_database_source_readback'
            or audit['database_receipt_sha256']!=sha(root/'receipt.json')):
        raise ValueError('Database readback mismatch')
    db=root/'domains.sqlite'
    if sha(db)!=receipt['artifacts']['domains.sqlite']:
        raise ValueError('Changed domain database')
    with Path(plan['nested_table']).open() as handle:
        nested={(r['containing_accession'],r['nested_accession']) for r in csv.DictReader(handle,delimiter='\t')}
    out=Path(plan['output'])
    if out.exists():
        raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free<plan['resources']['minimum_free_disk_gib']*2**30:
        raise RuntimeError('Insufficient disk')
    out.mkdir()
    start=time.time()
    conn=sqlite3.connect('file:'+str(db.resolve())+'?mode=ro',uri=True)
    conn.row_factory=sqlite3.Row
    fields=['hit_id','pfam_accession','pfam_clan','pfam_type','alignment_start','alignment_end',
            'envelope_start','envelope_end','hmm_coverage','independent_evalue','domain_score']
    sql=('SELECT q.sequence_id,q.search_partition,'+','.join('h.'+f for f in fields)
         +' FROM queries q LEFT JOIN hits h ON h.sequence_id=q.sequence_id ORDER BY q.sequence_id')
    totals=Counter()
    policies={coord+'_'+rank:Counter() for coord in ['alignment','envelope'] for rank in ['evalue','bitscore']}
    table=out/'query_competition.tsv.gz'
    def state(status):
        p=out/'state.tmp.json'
        p.write_text(json.dumps(dict(status=status,totals=dict(totals),elapsed_seconds=time.time()-start))+'\n')
        p.replace(out/'state.json')
    state('running')
    output_fields=['sequence_id','search_partition','raw_hits','raw_partial_hmm_hits_below_070',
                   'retained_sets_agree','status','policies_json']
    with gzip.open(table,'wt',compresslevel=3) as handle:
        writer=csv.DictWriter(handle,output_fields,delimiter='\t',lineterminator='\n')
        writer.writeheader()
        for sequence,group in itertools.groupby(conn.execute(sql),key=lambda r:r['sequence_id']):
            source=list(group)
            hits=[dict(r) for r in source if r['hit_id'] is not None]
            for h in hits:
                for field in ['alignment_start','alignment_end','envelope_start','envelope_end']:
                    h[field]=int(h[field])
            partial=sum(Decimal(h['hmm_coverage'])<Decimal('0.70') for h in hits)
            results={}
            for coord in ['alignment','envelope']:
                for rank in ['evalue','bitscore']:
                    key=coord+'_'+rank
                    result=compete(hits,nested,coord,rank)
                    if len(result['decisions'])!=len(hits):
                        raise ValueError('Incomplete hit disposition')
                    results[key]=result
                    pc=policies[key]
                    pc['retained_hits']+=len(result['retained_hits'])
                    pc['suppressed_hits']+=len(hits)-len(result['retained_hits'])
                    pc['queries_with_primary_rank_ties']+=bool(result['primary_rank_ties'])
                    pc['queries_with_unresolved_overlaps']+=bool(result['unresolved_overlap_pairs'])
                    pc['queries_with_candidate_nesting']+=bool(result['candidate_nested_pairs'])
            agree=len({tuple(sorted(r['retained_hits'])) for r in results.values()})==1
            status='no_GA_hit_not_proven_absence' if not hits else 'candidate_annotations_require_validation'
            writer.writerow(dict(sequence_id=sequence,search_partition=source[0]['search_partition'],
                                 raw_hits=len(hits),raw_partial_hmm_hits_below_070=partial,
                                 retained_sets_agree=int(agree),status=status,
                                 policies_json=json.dumps(results,separators=(',',':'),sort_keys=True)))
            totals.update(queries=1,hits=len(hits),partial_hmm_hits_below_070=partial,
                          queries_with_policy_disagreement=int(not agree),queries_without_hits=int(not hits))
            if totals['queries']%100000==0:
                handle.flush()
                state('running')
                print(json.dumps(dict(totals)),flush=True)
                if table.stat().st_size>plan['resources']['output_allowance_gib']*2**30:
                    raise RuntimeError('Output allowance exceeded')
                if shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30:
                    raise RuntimeError('Free disk gate reached')
    conn.close()
    if totals['queries']!=receipt['queries'] or totals['hits']!=receipt['total_hit_rows']:
        raise ValueError('Incomplete full-domain scope')
    for pc in policies.values():
        if pc['retained_hits']+pc['suppressed_hits']!=totals['hits']:
            raise ValueError('Incomplete policy dispositions')
    for path,digest in plan['pins'].items():
        if sha(path)!=digest:
            raise ValueError('Pin changed during execution')
    result=dict(status='completed_full_competition_requires_independent_readback',totals=dict(totals),
                policies={k:dict(v) for k,v in policies.items()},plan_sha256=sha(a.plan),
                script_sha256=sha(Path(__file__)),artifacts={table.name:sha(table)},
                elapsed_seconds=time.time()-start,
                scope='Every query and every raw hit has a disposition under all four policies. All raw annotation fields remain in the pinned source database. Retained hits are candidate annotations, not validated architectures or domain gains/losses. Independent empirical readback remains required.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    state(result['status'])
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    main()
