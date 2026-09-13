#!/usr/bin/env python3
"""Retrieve assembly-matched statistics with explicit metric scopes and provenance."""
import concurrent.futures
import csv
import datetime
import fcntl
import hashlib
import json
import math
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_statistics(text, accession):
    headers, values = {}, {}
    for line in text.splitlines():
        if line.startswith('# ') and ':' in line:
            key, value = line[2:].split(':', 1)
            headers[key.strip()] = value.strip()
        elif line.strip() and not line.startswith('#'):
            fields = line.split('\t')
            if len(fields) != 6:
                raise ValueError('Expected six assembly statistic fields')
            key = tuple(fields[:5])
            if key in values:
                raise ValueError('Duplicate assembly statistic')
            number = float(fields[5])
            if not math.isfinite(number) or number < 0:
                raise ValueError('Invalid assembly statistic value')
            values[key] = number
    accessions = {v.split()[0] for k, v in headers.items() if k in ['GenBank assembly accession', 'RefSeq assembly accession']}
    if accession not in accessions:
        raise ValueError('Statistics report belongs to another assembly version')
    all_stats = {k[4]: v for k, v in values.items() if k[:4] == ('all', 'all', 'all', 'all')}
    primary = {k[4]: v for k, v in values.items() if k[:4] == ('Primary Assembly', 'all', 'all', 'all')}
    if all_stats.get('total-length', 0) <= 0:
        raise ValueError('Missing positive whole-assembly length')
    for scope in [all_stats, primary]:
        if not scope:
            continue
        for field in ['ungapped-length', 'total-gap-length', 'contig-N50', 'scaffold-N50']:
            if field in scope and scope[field] > scope.get('total-length', float('inf')):
                raise ValueError('Component statistic exceeds assembly length')
        if 'gc-perc' in scope and scope['gc-perc'] > 100:
            raise ValueError('Invalid GC percentage')
    return {'headers': headers, 'all_assembly_metrics': all_stats, 'primary_assembly_metrics': primary,
            'metric_rows': len(values)}


def retrieve(row, folder):
    accession = row['assembly_accession']
    url = row['proteome_url'].removesuffix('_protein.faa.gz') + '_assembly_stats.txt'
    path = folder / (accession + '_assembly_stats.txt')
    receipt_path = folder / (accession + '.receipt.json')
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text())
        if receipt['assembly_accession'] != accession or receipt['url'] != url or sha(path) != receipt['sha256']:
            raise ValueError('Changed cached assembly report')
        parse_statistics(path.read_text(), accession)
        return receipt
    checksum_url = url.rsplit('/', 1)[0] + '/md5checksums.txt'
    for attempt in range(3):
        try:
            with urllib.request.urlopen(checksum_url, timeout=45) as response:
                checksum_bytes = response.read()
            name = url.rsplit('/', 1)[1]
            checksums = [line.split()[0] for line in checksum_bytes.decode().splitlines()
                         if len(line.split()) == 2 and line.split()[1].lstrip('./') == name]
            if len(checksums) != 1:
                raise ValueError('Missing or ambiguous publisher checksum')
            with urllib.request.urlopen(url, timeout=60) as response:
                data = response.read()
            if hashlib.md5(data).hexdigest() != checksums[0]:
                raise ValueError('Publisher checksum mismatch')
            parsed = parse_statistics(data.decode(), accession)
            break
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)
    temp = path.with_suffix('.partial')
    temp.write_bytes(data)
    temp.replace(path)
    result = {'taxon_id': row['taxon_id'], 'assembly_accession': accession, 'status': 'verified',
        'url': url, 'path': str(path.relative_to(ROOT)), 'sha256': sha(path), 'bytes': len(data),
        'publisher_md5': checksums[0], 'checksum_url': checksum_url,
        'checksum_listing_sha256': hashlib.sha256(checksum_bytes).hexdigest(),
        'retrieved_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), **parsed}
    temp = receipt_path.with_suffix('.partial')
    temp.write_text(json.dumps(result, indent=2) + '\n')
    temp.replace(receipt_path)
    return result


def main():
    manifest = ROOT / 'metadata/analysis_manifest.tsv'
    with manifest.open() as handle:
        taxa = list(csv.DictReader(handle, delimiter='\t'))
    folder = ROOT / 'data/assembly_statistics'
    folder.mkdir(parents=True, exist_ok=True)
    lock = (folder / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    ncbi = [r for r in taxa if r['proteome_url'].startswith('https://ftp.ncbi.nlm.nih.gov/')]
    records = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(retrieve, r, folder): r for r in ncbi}
        for future in concurrent.futures.as_completed(futures):
            row = futures[future]
            try:
                result = future.result()
            except Exception as error:
                result = {'taxon_id': row['taxon_id'], 'assembly_accession': row['assembly_accession'],
                          'status': 'failed', 'error': str(error)}
            records[row['taxon_id']] = result
            print(len(records), row['taxon_id'], result['status'], result.get('error', ''), flush=True)
    rows = []
    for taxon in taxa:
        r = records.get(taxon['taxon_id'], {'status': 'external_source_see_separate_metrics'})
        h = r.get('headers', {})
        row = {k: taxon[k] for k in ['taxon_id', 'species_name', 'study_role', 'assembly_accession']}
        row.update(statistics_status=r['status'], reported_organism=h.get('Organism name', ''),
                   reported_assembly_type=h.get('Assembly type', ''), genome_representation=h.get('Genome representation', ''),
                   assembly_method=h.get('Assembly method', ''), sequencing_technology=h.get('Sequencing technology', ''),
                   reported_coverage=h.get('Genome coverage', ''), biosample=h.get('BioSample', ''), bioproject=h.get('BioProject', ''))
        for label, field in [('all', 'all_assembly_metrics'), ('primary', 'primary_assembly_metrics')]:
            values = r.get(field, {})
            for metric in ['total-length', 'ungapped-length', 'total-gap-length', 'contig-count', 'contig-N50',
                           'contig-L50', 'scaffold-count', 'scaffold-N50', 'scaffold-L50', 'gc-perc']:
                row[label + '_' + metric] = values.get(metric, '')
            row[label + '_gap_fraction'] = values['total-gap-length'] / values['total-length'] if (
                'total-gap-length' in values and values.get('total-length', 0) > 0) else ''
        rows.append(row)
    table = ROOT / 'metadata/assembly_quality_metrics.tsv'
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    provenance = folder / 'retrieval_manifest.json'
    provenance.write_text(json.dumps([records[r['taxon_id']] for r in ncbi], indent=2) + '\n')
    result = {'manifest_sha256': sha(manifest), 'taxa': len(taxa), 'ncbi_requested': len(ncbi),
        'ncbi_verified': sum(r['status'] == 'verified' for r in records.values()),
        'ncbi_failed': sum(r['status'] == 'failed' for r in records.values()), 'external_sources_not_in_ncbi_reports': len(taxa) - len(ncbi),
        'downloaded_bytes': sum(r.get('bytes', 0) for r in records.values()),
        'provenance_path': str(provenance.relative_to(ROOT)), 'provenance_sha256': sha(provenance),
        'metrics_sha256': sha(table), 'script_sha256': sha(Path(__file__)),
        'interpretation': 'Deposited assembly metrics, not recomputed from FASTA; whole-assembly and primary-assembly scopes remain distinct. Missing fields are not zeros. Contiguity is not contamination assessment.'}
    (ROOT / 'metadata/assembly_quality_receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
