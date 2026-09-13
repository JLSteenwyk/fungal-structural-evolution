#!/usr/bin/env python3
"""Measure FASTA-record contiguity and base composition for all external genomes."""
import csv
import gzip
import json
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from retrieve_assembly_statistics import ROOT, sha


def n50(lengths):
    if not lengths or min(lengths) <= 0:
        raise ValueError('Positive nonempty sequence lengths required')
    total, cumulative = sum(lengths), 0
    for i, length in enumerate(sorted(lengths, reverse=True), 1):
        cumulative += length
        if cumulative * 2 >= total:
            return length, i
    raise AssertionError('Unreachable N50')


def main():
    manifest = ROOT / 'metadata/analysis_manifest.tsv'
    receipts_path = ROOT / 'metadata/external_genome_receipts.json'
    receipts = json.loads(receipts_path.read_text())
    with manifest.open() as handle:
        taxa = [r for r in csv.DictReader(handle, delimiter='\t')
                if not r['proteome_url'].startswith('https://ftp.ncbi.nlm.nih.gov/')]
    rows, sources = [], []
    for taxon in taxa:
        candidates = [r for r in receipts if r['url'] == taxon['genome_url']]
        if len(candidates) != 1:
            raise ValueError('Ambiguous external genome source')
        source = candidates[0]
        path = ROOT / source['path']
        if sha(path) != source['sha256']:
            raise ValueError('Changed external genome FASTA')
        counts, lengths, seen = Counter(), [], set()
        opener = gzip.open if path.suffix == '.gz' else open
        with opener(path, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                sequence = str(record.seq).upper()
                if record.id in seen or not sequence or set(sequence) - set('ACGTRYSWKMBDHVN'):
                    raise ValueError('Duplicate, empty or non-IUPAC DNA record')
                seen.add(record.id)
                lengths.append(len(sequence))
                counts.update(sequence)
        value, l50 = n50(lengths)
        canonical = sum(counts[b] for b in 'ACGT')
        if not canonical:
            raise ValueError('No canonical DNA bases')
        total = sum(lengths)
        row = {k: taxon[k] for k in ['taxon_id', 'species_name', 'study_role']}
        row.update(fasta_records=len(lengths), total_length=total, record_N50=value, record_L50=l50,
            longest_record=max(lengths), N_bases=counts['N'], N_fraction=counts['N'] / total,
            other_ambiguous_bases=total - canonical - counts['N'],
            gc_percent_canonical=100 * (counts['G'] + counts['C']) / canonical,
            scope='All deposited FASTA records; no gap splitting or organelle filtering', status='measured')
        rows.append(row)
        sources.append({'taxon_id': taxon['taxon_id'], 'path': source['path'], 'url': source['url'],
                        'sha256': source['sha256'], 'citation': source['citation'], 'article_version': source['article_version']})
        print(taxon['taxon_id'], len(lengths), value, flush=True)
    table = ROOT / 'metadata/external_assembly_quality.tsv'
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    receipt = {'taxa': len(rows), 'manifest_sha256': sha(manifest), 'source_receipts_sha256': sha(receipts_path),
        'sources': sources, 'table_sha256': sha(table), 'script_sha256': sha(Path(__file__)),
        'interpretation': 'Direct FASTA-record statistics; record N50 is not assumed equivalent to NCBI contig N50. Ambiguous bases and all deposited sequences remain explicit. No contamination assessment.'}
    (ROOT / 'metadata/external_assembly_quality_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
