#!/usr/bin/env python3
"""Freeze the live retrieval-log prefix and check all known marker gaps for cached models."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import shutil
from datetime import datetime, timezone


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Plan changed')
        for path,digest in plan['pins'].items():
            if sha(path)!=digest:raise ValueError('Pinned input changed: '+path)
    verify();gap=Path(plan['gaps']);gr=json.loads((gap/'receipt.json').read_text())
    for name in ['gap_unique_sequences.tsv','gap_marker_records.tsv']:
        if sha(gap/name)!=gr['artifacts'][name]:raise ValueError('Unbound gap table')
    with (gap/'gap_unique_sequences.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    missing={r['sequence_sha256']:r for r in rows}
    if len(missing)!=len(rows) or len(rows)!=gr['missing_unique_sequences']:raise ValueError('Gap scope differs')
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    if shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk')
    source=Path(plan['inventory']);frozen=out/'retrieval_prefix.jsonl'
    # Copy only the byte extent observed on opening. A concurrent partial final
    # line is excluded; all complete prior records are retained unchanged.
    with source.open('rb') as src, frozen.open('wb') as dst:
        size=source.stat().st_size;remaining=size
        while remaining:
            b=src.read(min(8*1024*1024,remaining))
            if not b:raise ValueError('Source truncated while freezing')
            dst.write(b);remaining-=len(b)
    with frozen.open('rb+') as f:
        end=f.seek(0,2)
        while end:
            f.seek(end-1)
            if f.read(1)==b'\n':break
            end-=1
        f.truncate(end)
    frozen_sha=sha(frozen);live_prefix=hashlib.sha256()
    with source.open('rb') as f:
        remaining=end
        while remaining:
            b=f.read(min(8*1024*1024,remaining))
            if not b:raise ValueError('Source prefix truncated')
            live_prefix.update(b);remaining-=len(b)
    if live_prefix.hexdigest()!=frozen_sha:raise ValueError('Live prefix changed during snapshot')
    latest={};records=0
    with frozen.open() as f:
        for line in f:
            record=json.loads(line);records+=1;accession=record['uniprot_accession']
            models=[m for m in record.get('models',[]) if m['sequence_sha256'] in missing] if record['status']=='verified' else []
            if models:latest[accession]=models
            else:latest.pop(accession,None)
    models={};accessions=defaultdict(set)
    for accession,entries in latest.items():
        for m in entries:
            key=m['model_id'],str(m['version'])
            if key in models and models[key]!=m:raise ValueError('Conflicting candidate records')
            models[key]=m;accessions[key].add(accession)
    counts=defaultdict(int);candidates=[]
    for key,m in sorted(models.items()):
        if int(m['length'])!=int(missing[m['sequence_sha256']]['length']) or sha(m['path'])!=m['sha256']:
            raise ValueError('Candidate sequence length or coordinate bytes differ')
        counts[m['sequence_sha256']]+=1;candidates.append(dict(model=m,uniprot_accessions=sorted(accessions[key])))
    target=out/'candidate_models.jsonl'
    with target.open('w') as f:
        for r in candidates:f.write(json.dumps(r,sort_keys=True)+'\n')
    with (out/'gap_sequence_availability.tsv').open('w') as f:
        fields=['sequence_sha256','length','verified_cached_candidates']
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for r in sorted(rows,key=lambda r:r['sequence_sha256']):
            w.writerow(dict(sequence_sha256=r['sequence_sha256'],length=r['length'],verified_cached_candidates=counts[r['sequence_sha256']]))
    with (gap/'gap_marker_records.tsv').open() as f:links=list(csv.DictReader(f,delimiter='\t'))
    if len(links)!=gr['missing_marker_records'] or {r['sequence_sha256'] for r in links}!=set(missing):raise ValueError('Gap link scope differs')
    verify()
    result=dict(status='complete_frozen_current_cache_marker_gap_screen',observed_at=datetime.now(timezone.utc).isoformat(),plan_sha256=ph,
                inventory_records=records,observed_live_bytes=size,frozen_complete_record_bytes=end,
                gap_unique_sequences=len(missing),gap_marker_records=len(links),candidate_models=len(models),
                gap_sequences_with_cached_candidates=sum(counts[k]>0 for k in missing),
                gap_marker_records_with_cached_candidates=sum(counts[r['sequence_sha256']]>0 for r in links),
                artifacts={p.name:sha(p) for p in out.iterdir()},
                scope='All prior frozen marker gaps checked against latest records per accession in a new immutable complete-line retrieval-log prefix. Candidate full-sequence hash/length and cached coordinate bytes checked. Not confidence qualification, a changed selected catalog, proof of public-database absence, or authorization to resume GPU prediction.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
