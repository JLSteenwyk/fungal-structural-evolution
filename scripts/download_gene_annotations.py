#!/usr/bin/env python3
"""Retrieve assembly-matched NCBI GFFs for gene/isoform reconciliation."""
import csv
import fcntl
import gzip
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def download(row):
    url = row['proteome_url'].removesuffix('_protein.faa.gz') + '_genomic.gff.gz'
    path = ROOT / 'data/annotations' / (row['assembly_accession'] + '_genomic.gff.gz')
    temp = path.with_suffix('.partial')
    result = {'taxon_id': row['taxon_id'], 'assembly_accession': row['assembly_accession'],
              'url': url, 'checked_at_utc': datetime.now(timezone.utc).isoformat()}
    try:
        checksum_url = url.rsplit('/', 1)[0] + '/md5checksums.txt'
        with urlopen(checksum_url, timeout=30) as response:
            checksums = response.read().decode()
        filename = url.rsplit('/', 1)[1]
        matches = [line.split()[0] for line in checksums.splitlines()
                   if len(line.split()) == 2 and line.split()[1].lstrip('./') == filename]
        if len(matches) != 1:
            raise ValueError('Missing or ambiguous publisher checksum')
        expected = matches[0]
        if not path.exists():
            with urlopen(url, timeout=90) as response, temp.open('wb') as out:
                for block in iter(lambda: response.read(1024 * 1024), b''):
                    out.write(block)
            if digest(temp, 'md5') != expected:
                raise ValueError('Downloaded annotation MD5 mismatch')
            temp.replace(path)
        if digest(path, 'md5') != expected:
            raise ValueError('Cached annotation MD5 mismatch')
        feature_counts = {}
        with gzip.open(path, 'rt') as handle:
            for line in handle:
                if line.startswith('##FASTA'):
                    break
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 9 or int(fields[3]) > int(fields[4]) or int(fields[3]) < 1:
                    raise ValueError('Invalid GFF feature row')
                feature_counts[fields[2]] = feature_counts.get(fields[2], 0) + 1
        if not feature_counts.get('CDS'):
            raise ValueError('Annotation contains no CDS features')
        result.update(status='validated', path=str(path.relative_to(ROOT)),
                      sha256=digest(path), publisher_md5=expected, checksum_url=checksum_url,
                      checksum_listing_sha256=hashlib.sha256(checksums.encode()).hexdigest(),
                      compressed_bytes=path.stat().st_size, feature_counts=feature_counts)
    except Exception as error:
        result.update(status='error', error=str(error))
    finally:
        temp.unlink(missing_ok=True)
    return result


def main():
    folder = ROOT / 'data/annotations'
    folder.mkdir(parents=True, exist_ok=True)
    lock = (folder / '.download.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    with (ROOT / 'metadata/analysis_manifest.tsv').open() as handle:
        rows = [r for r in csv.DictReader(handle, delimiter='\t')
                if r['proteome_url'].startswith('https://ftp.ncbi.nlm.nih.gov/')
                and r['proteome_url'].endswith('_protein.faa.gz')]
    cache = ROOT / 'data/raw/annotation_downloads.jsonl'
    done = {}
    if cache.exists():
        for line in cache.read_text().splitlines():
            record = json.loads(line)
            if record['status'] == 'validated':
                path = ROOT / record['path']
                if path.exists() and digest(path) == record['sha256']:
                    done[record['taxon_id']] = record
    pending = [r for r in rows if r['taxon_id'] not in done]
    print('Annotation candidates:', len(rows), 'pending:', len(pending), flush=True)
    with cache.open('a') as out, ThreadPoolExecutor(max_workers=2) as pool:
        for i, future in enumerate(as_completed([pool.submit(download, r) for r in pending]), 1):
            result = future.result()
            done[result['taxon_id']] = result
            out.write(json.dumps(result) + '\n')
            out.flush()
            print(i, result['taxon_id'], result['status'], result.get('error', ''), flush=True)
    (ROOT / 'metadata/annotation_download_receipts.json').write_text(
        json.dumps([done[r['taxon_id']] for r in rows], indent=2) + '\n')


if __name__ == '__main__':
    main()
