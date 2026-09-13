#!/usr/bin/env python3
"""Prepare all representative proteins for domain search, reusing exact marker sequences."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from prepare_pfam import ROOT, digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = ROOT / 'results/gene_representatives/full-v2'
    receipt = json.loads((source / 'receipt.json').read_text())
    manifest = ROOT / 'metadata/analysis_manifest.tsv'
    with manifest.open() as handle:
        expected = {r['taxon_id'] for r in csv.DictReader(handle, delimiter='\t')}
    if receipt['status'] != 'complete' or {r['taxon_id'] for r in receipt['taxa']} != expected:
        raise ValueError('Complete representative proteomes for the full taxon set required')
    marker_input = ROOT / 'data/domains/marker-inputs-v1'
    mr = json.loads((marker_input / 'receipt.json').read_text())
    if digest(marker_input / 'sequences.faa') != mr['artifacts']['sequences.faa']:
        raise ValueError('Changed marker-domain sequences')
    marker_sequences = {r.id[1:]: str(r.seq) for r in SeqIO.parse(marker_input / 'sequences.faa', 'fasta')}
    if any(hashlib.sha256(s.encode()).hexdigest() != h for h, s in marker_sequences.items()):
        raise ValueError('Marker sequence identity mismatch')
    if args.output.exists():
        raise FileExistsError('Use a new immutable full-domain input directory')
    args.output.mkdir(parents=True)
    unique, rows, reused, count, residues = {}, [], set(), 0, 0
    targets = args.output / 'additional_sequences.faa'
    links = args.output / 'protein_links.tsv'
    with targets.open('w') as out, links.open('w') as link_handle:
        writer = csv.writer(link_handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['taxon_id', 'protein_id', 'sequence_id', 'query_source'])
        for item in sorted(receipt['taxa'], key=lambda r: r['taxon_id']):
            path = ROOT / item['path']
            if digest(path) != item['sha256']:
                raise ValueError('Changed representative proteome')
            seen, n, marker_links, additional_links = set(), 0, 0, 0
            for record in SeqIO.parse(path, 'fasta'):
                if record.id in seen or not record.seq:
                    raise ValueError('Duplicate or empty representative protein')
                seen.add(record.id)
                sequence = str(record.seq)
                h = hashlib.sha256(sequence.encode()).hexdigest()
                if h in marker_sequences:
                    if marker_sequences[h] != sequence:
                        raise ValueError('Sequence hash collision')
                    partition = 'marker-inputs-v1'
                    reused.add(h)
                    marker_links += 1
                else:
                    if h in unique and unique[h] != sequence:
                        raise ValueError('Sequence hash collision')
                    if h not in unique:
                        unique[h] = sequence
                        residues += len(sequence)
                        out.write(f'>S{h}\n{sequence}\n')
                    partition = 'additional_full_proteome'
                    additional_links += 1
                writer.writerow([item['taxon_id'], record.id, 'S' + h, partition])
                n += 1
            if n != item['selected_proteins']:
                raise ValueError('Representative protein count differs from receipt')
            count += n
            rows.append({'taxon_id': item['taxon_id'], 'proteins': n, 'source_sha256': item['sha256'],
                         'marker_query_links': marker_links, 'additional_query_links': additional_links})
            print(len(rows), item['taxon_id'], n, flush=True)
    result = {'status': 'complete_search_input_preparation', 'taxa': len(rows), 'representative_proteins': count,
        'additional_unique_sequences': len(unique), 'additional_unique_residues': residues,
        'reused_unique_marker_sequences': len(reused), 'full_proteome_unique_sequences': len(unique) + len(reused),
        'marker_query_protein_links': sum(r['marker_query_links'] for r in rows),
        'source_receipt_sha256': digest(source / 'receipt.json'), 'manifest_sha256': digest(manifest),
        'marker_input_receipt_sha256': digest(marker_input / 'receipt.json'), 'taxon_inputs': rows,
        'script_sha256': digest(Path(__file__)), 'artifacts': {p.name: digest(p) for p in args.output.iterdir()},
        'interpretation': 'Exact sequence deduplication for domain search only. Every selected protein/taxon link retained, including unresolved gene mappings. Marker results may be reused only after that complete search passes validation. Additional full-proteome domain searches have not been executed by this script.'}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'taxon_inputs'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
