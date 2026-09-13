#!/usr/bin/env python3
"""Acquire publisher-verified FCS reports for exact selected assembly versions."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import time
import urllib.request
import urllib.error

FIELDS=['seq_id','start_pos','end_pos','seq_len','action','contam_type','coverage','contam_details']


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def parse_report(text):
    lines=text.splitlines()
    meta=[json.loads(s[2:]) for s in lines if s.startswith('##')]
    headers=[s[1:].split('\t') for s in lines if s.startswith('#') and not s.startswith('##')]
    if len(meta)!=1 or meta[0][0][0]!='FCS genome report' or headers!=[FIELDS]: raise ValueError('Unrecognized FCS header')
    rows=[];seen=set();lengths={}
    for line in lines:
        if not line.strip() or line.startswith('#'):continue
        parts=line.split('\t')
        if len(parts)!=8:raise ValueError('Wrong FCS field count')
        row=dict(zip(FIELDS,parts)); key=tuple(parts)
        if key in seen:raise ValueError('Duplicate FCS record')
        seen.add(key)
        for name in ['start_pos','end_pos','seq_len']:row[name]=int(row[name])
        if not 1<=row['start_pos']<=row['end_pos']<=row['seq_len']:raise ValueError('Invalid one-based inclusive interval')
        cov=float(row['coverage'])
        if not math.isfinite(cov) or not 0<=cov<=100:raise ValueError('Invalid reported coverage')
        if row['seq_id'] in lengths and lengths[row['seq_id']]!=row['seq_len']:raise ValueError('Conflicting sequence lengths')
        lengths[row['seq_id']]=row['seq_len'];rows.append(row)
    return meta[0],rows


def union_bases(rows, actions):
    grouped=defaultdict(list)
    for r in rows:
        if r['action'] in actions:grouped[r['seq_id']].append((r['start_pos'],r['end_pos']))
    total=0
    for intervals in grouped.values():
        end=0
        for start,stop in sorted(intervals):
            total+=max(0,stop-max(end,start-1));end=max(end,stop)
    return total


def fetch(url):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url,timeout=40) as response:return response.read()
        except urllib.error.HTTPError as exc:
            if exc.code==404 or attempt==2:raise
        except OSError:
            if attempt==2:raise
        time.sleep(2**attempt)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use an immutable output directory')
    with a.manifest.open() as f:taxa=list(csv.DictReader(f,delimiter='\t'))
    if len(taxa)!=len({r['taxon_id'] for r in taxa}):raise ValueError('Duplicate manifest taxon')
    a.output.mkdir(parents=True);(a.output/'reports').mkdir();manifest_hash=sha(a.manifest)
    def retrieve(row):
        accession=row['assembly_accession'];base={'taxon_id':row['taxon_id'],'assembly_accession':accession,'retrieved_utc':datetime.now(timezone.utc).isoformat()}
        if not row['proteome_url'].startswith('https://ftp.ncbi.nlm.nih.gov/'):
            return dict(base,status='external_source_no_ncbi_report_requested'),[],{}
        url=row['proteome_url'].removesuffix('_protein.faa.gz')+'_fcs_report.txt';name=url.rsplit('/',1)[1]
        if not name.startswith(accession+'_'):raise ValueError('URL assembly version differs')
        folder=a.output/'reports'/accession;folder.mkdir();base.update(url=url,checksum_url=url.rsplit('/',1)[0]+'/md5checksums.txt')
        try:
            listing=fetch(base['checksum_url']);(folder/'md5checksums.txt').write_bytes(listing)
            hashes=[s.split()[0] for s in listing.decode().splitlines() if len(s.split())==2 and s.split()[1].removeprefix('./')==name]
            raw=fetch(url);(folder/name).write_bytes(raw)
            if len(hashes)!=1 or hashlib.md5(raw).hexdigest()!=hashes[0]:raise ValueError('Missing or mismatched publisher checksum')
            meta,records=parse_report(raw.decode());base.update(status='publisher_verified_report',publisher_md5=hashes[0],report_sha256=sha(folder/name),checksum_listing_sha256=sha(folder/'md5checksums.txt'),report_path=str(folder/name),report_rows=len(records),header=meta)
            joined=[dict(taxon_id=row['taxon_id'],assembly_accession=accession,**r) for r in records]
            info={'reported_actions':dict(Counter(r['action'] for r in records)),'flagged_sequence_count':len({r['seq_id'] for r in records}),'exclude_fix_trim_union_bp':union_bases(records,{'EXCLUDE','FIX','TRIM'}),'review_rare_union_bp':union_bases(records,{'REVIEW_RARE'}),'review_union_bp':union_bases(records,{'REVIEW'}),'unknown_actions':sorted({r['action'] for r in records}-{'EXCLUDE','FIX','TRIM','REVIEW_RARE','REVIEW','INFO','MITOCHONDRION','PLASTID'})}
        except Exception as exc:
            base.update(status='report_http_404' if isinstance(exc,urllib.error.HTTPError) and exc.code==404 else 'report_unresolved_error',error_type=type(exc).__name__,error=str(exc));joined=[];info={}
        base['saved_artifacts']={p.name:sha(p) for p in folder.iterdir() if p.is_file()}
        (folder/'receipt.json').write_text(json.dumps(base,indent=2)+'\n')
        return base,joined,info
    receipts=[];all_rows=[];summaries=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures={pool.submit(retrieve,row):row for row in taxa}
        for i,future in enumerate(as_completed(futures),1):
            source=futures[future];r,rows,info=future.result();receipts.append(r);all_rows.extend(rows)
            summary={'taxon_id':source['taxon_id'],'assembly_accession':source['assembly_accession'],'status':r['status'],'report_rows':r.get('report_rows',''),'flagged_sequence_count':info.get('flagged_sequence_count',''),'exclude_fix_trim_union_bp':info.get('exclude_fix_trim_union_bp',''),'review_rare_union_bp':info.get('review_rare_union_bp',''),'review_union_bp':info.get('review_union_bp',''),'reported_actions':json.dumps(info.get('reported_actions',{}),sort_keys=True),'unknown_actions':';'.join(info.get('unknown_actions',[])),'report_path':r.get('report_path',''),'report_sha256':r.get('report_sha256','')}
            summaries.append(summary)
            if i%50==0:print('Retrieved/reviewed',i,'of',len(taxa),flush=True)
    for name,rows,fields in [('taxon_summary.tsv',summaries,list(summaries[0])),('flagged_regions.tsv',all_rows,['taxon_id','assembly_accession']+FIELDS)]:
        rows.sort(key=lambda r:(r['taxon_id'],str(r.get('seq_id','')),int(r.get('start_pos',0))))
        with (a.output/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    if sha(a.manifest)!=manifest_hash:raise ValueError('Manifest changed during acquisition')
    receipt={'status':'complete_selected_assembly_fcs_report_inventory','manifest_sha256':manifest_hash,'script_sha256':sha(Path(__file__)),'taxa':len(taxa),'status_counts':dict(Counter(r['status'] for r in summaries)),'verified_reports_with_rows':sum(r['status']=='publisher_verified_report' and r['report_rows']>0 for r in summaries),'flagged_region_rows':len(all_rows),'taxa_with_exclude_fix_trim_regions':sum(isinstance(r['exclude_fix_trim_union_bp'],int) and r['exclude_fix_trim_union_bp']>0 for r in summaries),'source_receipts':sorted(receipts,key=lambda r:r['taxon_id']),'artifacts':{n:sha(a.output/n) for n in ['taxon_summary.tsv','flagged_regions.tsv']},'interpretation':'Publisher FCS observations for exact assembly-version URLs, not independent rescreening or proof of contamination absence. EXCLUDE/FIX/TRIM, REVIEW_RARE, REVIEW, INFO and organelle records remain distinct. Missing/error/external reports are not zero contamination. Gene/CDS overlap and biological review remain pending; no sequences removed.'}
    (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_receipts','artifacts']},indent=2))


if __name__=='__main__':main()
