#!/usr/bin/env python3
"""Acquire exact assembly DNA for every selected taxon with publisher-bound provenance."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import threading
import time
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from reference_measurement_union_sources import bind, verify


def digest(path, algorithm='sha256'):
    value = hashlib.new(algorithm)
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b''): value.update(chunk)
    return value.hexdigest()


def publisher_md5(listing, filename):
    found = []
    for line in listing.decode('utf-8').splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[1].lstrip('./') == filename:
            assert re.fullmatch('[0-9a-fA-F]{32}', fields[0])
            found.append(fields[0].lower())
    if len(found) != 1: raise ValueError('Missing or ambiguous exact publisher genome checksum')
    return found[0]


def read_fasta(path, contigs):
    """Stream original DNA, retaining case/ambiguity and every deposited record."""
    counts, lengths, seen = Counter(), [], set()
    identifier = None; length = 0; lower = 0; whitespace = 0
    original_hash = None; canonical_hash = None
    alphabet = set(b'ACGTRYSWKMBDHVN')
    opener = gzip.open if Path(path).suffix == '.gz' else open
    with opener(path, 'rb') as handle, Path(contigs).open('x') as out:
        def finish():
            if identifier is None: return
            if length == 0: raise ValueError('Empty genomic FASTA record')
            lengths.append(length)
            out.write(json.dumps(dict(sequence_id=identifier, description=description, length=length,
                                      sequence_sha256=original_hash.hexdigest(),
                                      uppercase_sequence_sha256=canonical_hash.hexdigest()),
                                 separators=(',', ':'), allow_nan=False)+'\n')
        for raw in handle:
            line = raw.rstrip(b'\r\n')
            if not line.strip(): continue
            if line.startswith(b'>'):
                finish(); description = line[1:].decode('utf-8')
                fields = description.split()
                if not fields or fields[0] in seen: raise ValueError('Empty or duplicate genomic FASTA identifier')
                identifier = fields[0]; seen.add(identifier); length = 0
                original_hash = hashlib.sha256(); canonical_hash = hashlib.sha256()
            else:
                if identifier is None: raise ValueError('DNA before genomic FASTA header')
                sequence = b''.join(line.split()); whitespace += len(line)-len(sequence)
                upper = sequence.upper()
                if set(upper)-alphabet: raise ValueError('Non-IUPAC genomic DNA')
                lower += len(sequence.translate(None, b'ACGTRYSWKMBDHVN'))
                length += len(sequence); counts.update(upper)
                original_hash.update(sequence); canonical_hash.update(upper)
        finish()  # Reading through EOF validates the gzip CRC/trailer too.
    if not lengths: raise ValueError('No genomic FASTA records')
    canonical = sum(counts[base] for base in b'ACGT')
    if not canonical: raise ValueError('No canonical genomic DNA bases')
    total = sum(lengths); accumulated = 0
    for rank, size in enumerate(sorted(lengths, reverse=True), 1):
        accumulated += size
        if 2*accumulated >= total: n50, l50 = size, rank; break
    return dict(fasta_records=len(lengths), total_length=total, longest_record=max(lengths),
                record_N50=n50, record_L50=l50, N_bases=counts[ord('N')],
                N_fraction=counts[ord('N')]/total, other_ambiguous_bases=total-canonical-counts[ord('N')],
                gc_percent_canonical=100*(counts[ord('G')]+counts[ord('C')])/canonical,
                lowercase_bases=lower, removed_sequence_whitespace=whitespace,
                base_counts={chr(k):v for k,v in sorted(counts.items())},
                scope='All deposited FASTA records; no gap splitting, organelle exclusion, taxonomic admission or contamination inference')


class Budget:
    def __init__(self, root, maximum_bytes, reserve_bytes):
        self.root = root; self.maximum = maximum_bytes; self.reserve = reserve_bytes
        self.bytes = 0; self.lock = threading.Lock()
    def account(self, size):
        with self.lock:
            if shutil.disk_usage(self.root).free < self.reserve:
                raise RuntimeError('Emergency genome disk reserve reached')
            if self.bytes+size > self.maximum:
                raise RuntimeError('Declared genome download-byte allowance reached')
            self.bytes += size


def acquire(entry, output, budget, per_file_bytes):
    taxon = entry['taxon']; identifier = taxon['taxon_id']
    assert re.fullmatch('[A-Za-z0-9_.-]+', identifier)
    root = Path(output)/identifier; root.mkdir(exist_ok=False)
    report = dict(taxon=taxon, source=entry['source'], attempted_http_requests=0)
    if entry['source'] == 'existing_external_publisher_bundle':
        original = entry['external_receipt']; path = Path(original['path']).resolve()
        assert digest(path) == original['sha256'] and digest(path, 'md5') == original['publisher_md5']
        report.update(path=str(path), sha256=original['sha256'], publisher_md5=original['publisher_md5'],
                      url=original['url'], bytes=path.stat().st_size, external_receipt=original,
                      new_download_bytes=0)
    else:
        assert entry['source'] == 'ncbi_exact_assembly_version'
        url = taxon['genome_url']; accession = taxon['assembly_accession']
        assert urlparse(url).hostname == 'ftp.ncbi.nlm.nih.gov'
        assert url == taxon['proteome_url'].removesuffix('_protein.faa.gz')+'_genomic.fna.gz'
        filename = url.rsplit('/', 1)[1]
        assert Path(filename).name == filename and filename.startswith(accession+'_')
        checksum_url = url.rsplit('/', 1)[0]+'/md5checksums.txt'
        downloaded = 0
        for attempt in range(1, 4):
            candidate = root/f'attempt_{attempt}.fna.gz'; listing_path = root/f'attempt_{attempt}.md5checksums.txt'
            try:
                request = Request(checksum_url, headers={'User-Agent':'fungal-structural-evolution/1.0'})
                report['attempted_http_requests'] += 1
                with urlopen(request, timeout=60) as response: listing = response.read(16*1024*1024+1)
                if len(listing) > 16*1024*1024: raise ValueError('Publisher checksum listing exceeds bound')
                listing_path.write_bytes(listing); expected = publisher_md5(listing, filename)
                request = Request(url, headers={'User-Agent':'fungal-structural-evolution/1.0'})
                report['attempted_http_requests'] += 1
                with urlopen(request, timeout=120) as response, candidate.open('xb') as handle:
                    content_length = response.headers.get('Content-Length'); length = 0
                    for chunk in iter(lambda: response.read(1024*1024), b''):
                        if length+len(chunk) > per_file_bytes: raise RuntimeError('Per-genome byte allowance reached')
                        budget.account(len(chunk)); handle.write(chunk); length += len(chunk); downloaded += len(chunk)
                if content_length is not None and int(content_length) != length:
                    raise ValueError('Genome response length mismatch')
                if digest(candidate, 'md5') != expected: raise ValueError('Publisher genomic FASTA checksum mismatch')
                path = candidate; break
            except RuntimeError: raise  # Budget/resource failures never become individual biological exclusions.
            except Exception:
                if attempt == 3: raise
                time.sleep(2**(attempt-1))
        report.update(path=str(path), sha256=digest(path), publisher_md5=expected, url=url,
                      bytes=path.stat().st_size, checksum_url=checksum_url,
                      checksum_listing=str(listing_path), checksum_listing_sha256=digest(listing_path),
                      http_content_length=content_length, new_download_bytes=downloaded, http_attempt=attempt)
    contigs = root/'contig_metadata.jsonl'; stats = read_fasta(path, contigs)
    report.update(status='verified_publisher_bound_genomic_dna', fasta=stats,
                  contig_metadata=str(contigs), contig_metadata_sha256=digest(contigs),
                  checked_utc=datetime.now(timezone.utc).isoformat())
    if entry['source'] == 'ncbi_exact_assembly_version':
        deposited = entry['statistics_receipt']; assert deposited['assembly_accession'] == taxon['assembly_accession']
        total = deposited['all_assembly_metrics']['total-length']; assert total == int(total)
        report['deposited_whole_assembly_length'] = int(total)
        report['fasta_minus_deposited_whole_assembly_length'] = stats['total_length']-int(total)
        report['assembly_scope_length_status'] = ('matches_deposited_all_scope' if stats['total_length'] == int(total)
                                                  else 'requires_assembly_scope_review')
    target = root/'receipt.json'
    with target.open('x') as handle: json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    return dict(**report, receipt_path=str(target), receipt_sha256=digest(target))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True); parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text()); pins = dict(plan['pins']); verify(pins)
    root = Path(plan['output']).resolve(); assert not root.exists()
    assert shutil.disk_usage(root.parent).free >= plan['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(); assert len(plan['entries']) == 526
    assert len({entry['taxon']['taxon_id'] for entry in plan['entries']}) == 526
    assert Counter(entry['source'] for entry in plan['entries']) == {'ncbi_exact_assembly_version':519,'existing_external_publisher_bundle':7}
    budget = Budget(root, plan['resources']['maximum_download_bytes'], plan['resources']['emergency_free_disk_gib']*2**30)
    records, counts, length_differences = [], Counter(), Counter(); started = time.monotonic()
    with (root/'taxon_dispositions.jsonl').open('x') as out, ThreadPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(acquire, entry, root, budget, plan['resources']['per_genome_bytes']): entry for entry in plan['entries']}
        for future in as_completed(futures):
            entry = futures[future]
            try: report = future.result()
            except RuntimeError:
                for pending in futures: pending.cancel()
                raise
            except Exception as error:
                report = dict(taxon=entry['taxon'], source=entry['source'], status='genome_acquisition_or_format_error',
                              error_type=type(error).__name__, error=str(error))
            records.append(report); counts[report['status']] += 1
            if 'assembly_scope_length_status' in report: length_differences[report['assembly_scope_length_status']] += 1
            out.write(json.dumps(report, separators=(',', ':'), allow_nan=False)+'\n'); out.flush()
            tmp = root/'state.tmp'
            tmp.write_text(json.dumps(dict(completed_taxa=len(records),expected_taxa=526,counts=dict(counts),
                                          assembly_scope_length_dispositions=dict(length_differences),
                                          total_download_bytes=budget.bytes,http_workers=2,
                                          elapsed_seconds=time.monotonic()-started),indent=2)+'\n'); tmp.replace(root/'state.json')
    assert len(records) == 526 and {r['taxon']['taxon_id'] for r in records} == {e['taxon']['taxon_id'] for e in plan['entries']}
    for report in records:
        if report['status'] != 'verified_publisher_bound_genomic_dna': continue
        for path_key, hash_key in (('path','sha256'),('contig_metadata','contig_metadata_sha256'),
                                   ('receipt_path','receipt_sha256'),('checksum_listing','checksum_listing_sha256')):
            if path_key in report: bind(pins,report[path_key],report[hash_key])
    for path in (args.plan,root/'taxon_dispositions.jsonl',Path(__file__)): bind(pins,path,digest(path))
    verify(pins)
    result = dict(status='completed_all_selected_assembly_dna_attempts_pending_independent_readback',
                  checked_utc=datetime.now(timezone.utc).isoformat(), expected_taxa=526, counts=dict(counts),
                  assembly_scope_length_dispositions=dict(length_differences), total_download_bytes=budget.bytes,
                  verified_genome_bases=sum(r['fasta']['total_length'] for r in records if 'fasta' in r),
                  source_hashes=pins, scientific_eligibility=False,gpu=False,new_predictions=0,new_cost_usd=0,
                  scope='Every selected taxon attempted, exact recorded NCBI assembly DNA/checksum listing or '
                        'existing exact publisher bundle retained. Complete original FASTA records and case/ambiguity '
                        'with per-contig sequence hashes; no gaps split, sequences removed or taxon labels changed. '
                        'Whole deposited length discrepancies are review dispositions, not silent substitutes. '
                        'Completion of attempts does not mean universal availability, independent full readback, '
                        'annotation/GFF/CDS correctness, absence of haplotigs/contamination, duplication or biological acceptance.')
    with args.receipt.open('x') as handle: json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__ == '__main__': main()
