#!/usr/bin/env python3
"""Verify gap cache screening using reverse-log latest-record selection."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import mmap
from pathlib import Path


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',required=True,type=Path)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text());rh=sha(rp)
    if receipt['status']!='complete_frozen_current_cache_marker_gap_screen' or receipt['plan_sha256']!=sha(args.plan):raise ValueError('Wrong producer receipt')
    for path,digest in plan['pins'].items():
        if sha(path)!=digest:raise ValueError('Changed input')
    for name,digest in receipt['artifacts'].items():
        if sha(root/name)!=digest:raise ValueError('Changed result')
    with (Path(plan['gaps'])/'gap_unique_sequences.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    missing={r['sequence_sha256']:int(r['length']) for r in rows}
    if len(rows)!=len(missing):raise ValueError('Repeated source sequence')
    seen=set();models={};accessions=defaultdict(set);records=0
    with (root/'retrieval_prefix.jsonl').open('rb') as f, mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as data:
        end=len(data)
        if end!=receipt['frozen_complete_record_bytes'] or data[end-1:end]!=b'\n':raise ValueError('Snapshot extent differs')
        while end:
            start=data.rfind(b'\n',0,end-1)+1;r=json.loads(data[start:end]);end=start;records+=1
            accession=r['uniprot_accession']
            if accession in seen:continue
            seen.add(accession)
            if r['status']!='verified':continue
            for m in r.get('models',[]):
                if m['sequence_sha256'] not in missing:continue
                key=m['model_id'],str(m['version'])
                if key in models and models[key]!=m:raise ValueError('Conflicting source model')
                models[key]=m;accessions[key].add(accession)
    expected=[dict(model=models[k],uniprot_accessions=sorted(accessions[k])) for k in sorted(models)]
    with (root/'candidate_models.jsonl').open() as f:actual=[json.loads(line) for line in f]
    if actual!=expected:raise ValueError('Candidates differ from reverse replay')
    counts=Counter()
    for m in models.values():
        digest=m['sequence_sha256']
        if m['length']!=missing[digest] or sha(m['path'])!=m['sha256']:raise ValueError('Invalid candidate bytes or length')
        counts[digest]+=1
    with (root/'gap_sequence_availability.tsv').open() as f:actual=list(csv.DictReader(f,delimiter='\t'))
    expected=[dict(sequence_sha256=k,length=str(missing[k]),verified_cached_candidates=str(counts[k])) for k in sorted(missing)]
    if actual!=expected:raise ValueError('Incomplete or wrong sequence table')
    with (Path(plan['gaps'])/'gap_marker_records.tsv').open() as f:links=list(csv.DictReader(f,delimiter='\t'))
    totals=dict(inventory_records=records,gap_unique_sequences=len(missing),gap_marker_records=len(links),candidate_models=len(models),gap_sequences_with_cached_candidates=sum(counts[k]>0 for k in missing),gap_marker_records_with_cached_candidates=sum(counts[r['sequence_sha256']]>0 for r in links))
    if any(receipt[k]!=v for k,v in totals.items()) or sha(rp)!=rh:raise ValueError('Summary or receipt differs')
    for name,digest in receipt['artifacts'].items():
        if sha(root/name)!=digest:raise ValueError('Result changed during readback')
    result=dict(status='passed_full_reverse_log_marker_gap_cache_readback',**totals,producer_receipt_sha256=rh,script_sha256=sha(__file__),scope='Reverse replay of every frozen log record, first occurrence per accession; all candidate records and all gap sequence counts reconstructed. Same source log and JSON format; not independent database retrieval or structural confidence qualification.')
    (root/'readback.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
