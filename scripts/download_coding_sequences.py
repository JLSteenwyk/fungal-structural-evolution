#!/usr/bin/env python3
"""Retrieve and validate assembly-matched NCBI CDS FASTAs for the full design."""
import csv
import fcntl
import gzip
import hashlib
import json
import re
import signal
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen
from Bio import SeqIO
from download_gene_annotations import digest

ROOT = Path(__file__).resolve().parents[1]


def validate_fasta(path):
    identifiers = set()
    bases = missing_protein_id = non_triplet = ambiguous = 0
    with gzip.open(path, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id in identifiers: raise ValueError('Duplicate CDS record identifier')
            identifiers.add(record.id)
            sequence = str(record.seq).upper()
            if not sequence or not set(sequence) <= set('ACGTRYSWKMBDHVN'):
                raise ValueError('Empty or non-IUPAC-DNA CDS sequence')
            bases += len(sequence)
            non_triplet += len(sequence) % 3 != 0
            ambiguous += not set(sequence) <= set('ACGT')
            missing_protein_id += not bool(re.search(r'\[protein_id=[^\]]+\]',record.description))
    if not identifiers: raise ValueError('No CDS records')
    return {'cds_records':len(identifiers), 'nucleotide_bases':bases,
        'records_without_protein_id':missing_protein_id,'non_triplet_length_records':non_triplet,
        'records_with_ambiguous_bases':ambiguous}


def download(row):
    folder=ROOT/'data/cds';url=row['cds_url']
    target=folder/(row['assembly_accession']+'_cds_from_genomic.fna.gz')
    temp=target.with_suffix('.partial')
    listing=folder/(row['assembly_accession']+'_md5checksums.txt')
    result={'taxon_id':row['taxon_id'],'assembly_accession':row['assembly_accession'],
        'url':url,'checked_at_utc':datetime.now(timezone.utc).isoformat()}
    try:
        checksum_url=url.rsplit('/',1)[0]+'/md5checksums.txt'
        with urlopen(checksum_url,timeout=60) as response: raw=response.read()
        filename=url.rsplit('/',1)[1]
        matches=[line.split()[0] for line in raw.decode().splitlines() if len(line.split())==2 and line.split()[1].lstrip('./')==filename]
        if len(matches)!=1: raise ValueError('Missing or ambiguous CDS publisher checksum')
        expected=matches[0]
        if listing.exists() and listing.read_bytes()!=raw:
            raise ValueError('Publisher checksum listing changed; preserve previous snapshot for review')
        listing.write_bytes(raw)
        if not target.exists():
            with urlopen(url,timeout=120) as response,temp.open('wb') as handle:
                for block in iter(lambda:response.read(1024*1024),b''):handle.write(block)
            if digest(temp,'md5')!=expected:raise ValueError('Downloaded CDS MD5 mismatch')
            temp.replace(target)
        if digest(target,'md5')!=expected:raise ValueError('Cached CDS MD5 mismatch')
        statistics=validate_fasta(target)
        result.update(status='validated_cds_fasta',path=str(target.relative_to(ROOT)),sha256=digest(target),
            compressed_bytes=target.stat().st_size,publisher_md5=expected,checksum_url=checksum_url,
            checksum_listing_path=str(listing.relative_to(ROOT)),checksum_listing_sha256=digest(listing),**statistics)
    except Exception as error:result.update(status='error',error=str(error))
    finally:temp.unlink(missing_ok=True)
    return result


def main():
    folder=ROOT/'data/cds';folder.mkdir(parents=True,exist_ok=True)
    lock=(folder/'.download.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    manifest=ROOT/'metadata/analysis_manifest.tsv'
    with manifest.open() as handle: all_rows=list(csv.DictReader(handle,delimiter='\t'))
    rows=[r for r in all_rows if r['cds_url'].startswith('https://ftp.ncbi.nlm.nih.gov/')]
    if len({r['taxon_id'] for r in rows})!=len(rows) or len({r['assembly_accession'] for r in rows})!=len(rows):raise ValueError('Duplicate taxon or assembly')
    for row in rows:
        if row['cds_url']!=row['proteome_url'].removesuffix('_protein.faa.gz')+'_cds_from_genomic.fna.gz':raise ValueError('CDS/protein assembly source differs')
    cache=ROOT/'data/raw/cds_downloads.jsonl';latest={}
    if cache.exists():
        for line in cache.read_text().splitlines():
            r=json.loads(line);latest[r['taxon_id']]=r
    pending=[];done={}
    for row in rows:
        r=latest.get(row['taxon_id'])
        if r and r['status']=='validated_cds_fasta' and r['url']==row['cds_url'] and r['assembly_accession']==row['assembly_accession']:
            if digest(ROOT/r['path'])!=r['sha256'] or digest(ROOT/r['checksum_listing_path'])!=r['checksum_listing_sha256']:raise ValueError('Changed cached validated source')
            done[row['taxon_id']]=r
        else:pending.append(row)
    stop=[False]
    def request_stop(signum,frame):stop[0]=True
    signal.signal(signal.SIGTERM,request_stop);signal.signal(signal.SIGINT,request_stop)
    print('NCBI CDS candidates',len(rows),'pending',len(pending),'other-source taxa',len(all_rows)-len(rows),flush=True)
    iterator=iter(pending);completed=0
    with cache.open('a') as handle,ThreadPoolExecutor(max_workers=2) as pool:
        active={}
        for _ in range(2):
            row=next(iterator,None)
            if row is not None:active[pool.submit(download,row)]=row['taxon_id']
        while active:
            ready,_=wait(active,return_when=FIRST_COMPLETED)
            for future in ready:
                active.pop(future);result=future.result();done[result['taxon_id']]=result
                handle.write(json.dumps(result)+'\n');handle.flush();completed+=1
                print(completed,result['taxon_id'],result['status'],result.get('error',''),flush=True)
                if not stop[0]:
                    row=next(iterator,None)
                    if row is not None:active[pool.submit(download,row)]=row['taxon_id']
    summary={'status':'complete_ncbi_cds_acquisition' if len(done)==len(rows) and all(r['status']=='validated_cds_fasta' for r in done.values()) else 'incomplete_ncbi_cds_acquisition',
        'planned_ncbi_taxa':len(rows),'validated_taxa':sum(r['status']=='validated_cds_fasta' for r in done.values()),
        'pending_or_error_taxa':[r['taxon_id'] for r in rows if done.get(r['taxon_id'],{}).get('status')!='validated_cds_fasta'],
        'other_source_taxa':[r['taxon_id'] for r in all_rows if r not in rows],
        'manifest_sha256':digest(manifest),'script_sha256':digest(Path(__file__)), 'interrupted':stop[0],
        'interpretation':'Publisher MD5, local SHA256, gzip integrity and CDS FASTA identity/alphabet checked. Missing protein IDs, non-triplet lengths and ambiguity are recorded, not discarded. This is not translation validation or codon-alignment/selection eligibility. Seven external-source taxa require their separately recorded extraction/acquisition paths.',
        'taxa':[done[r['taxon_id']] for r in rows if r['taxon_id'] in done]}
    (ROOT/'metadata/cds_download_receipts.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary['status'],summary['validated_taxa'],flush=True)


if __name__=='__main__':main()
