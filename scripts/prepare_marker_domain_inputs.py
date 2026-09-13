#!/usr/bin/env python3
"""Prepare exact-sequence-deduplicated domain inputs with all marker/taxon links."""
import csv
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from prepare_pfam import ROOT, digest


def main():
    source = ROOT / 'results/phylogeny/markers-full-v1'
    receipt = json.loads((source / 'receipt.json').read_text())
    if receipt['status'] != 'complete_extraction' or receipt['pending_taxa']:
        raise ValueError('Complete marker extraction required')
    output = ROOT / 'data/domains/marker-inputs-v1'
    if output.exists():
        raise FileExistsError('Immutable marker domain inputs already exist')
    with (source / 'protein_mapping.tsv').open() as handle:
        mapping = list(csv.DictReader(handle, delimiter='\t'))
    expected = {(r['marker'], r['taxon_id']): r for r in mapping}
    if len(expected) != len(mapping):
        raise ValueError('Duplicate marker/taxon source mapping')
    with (source / 'occupancy.tsv').open() as handle:
        files = list(csv.DictReader(handle, delimiter='\t'))
    unique, seen = {}, set()
    for row in files:
        path = source / 'unaligned' / (row['marker'] + '.faa')
        if digest(path) != row['sha256']:
            raise ValueError('Changed extracted marker sequence file')
        for record in SeqIO.parse(path, 'fasta'):
            key = (row['marker'], record.id)
            sequence = str(record.seq).upper()
            h = hashlib.sha256(sequence.encode()).hexdigest()
            if key in seen or h != expected[key]['sequence_sha256']:
                raise ValueError('Duplicate or changed marker sequence')
            if h in unique and unique[h] != sequence:
                raise ValueError('Sequence hash collision')
            unique[h] = sequence
            seen.add(key)
    if seen != expected.keys() or len(seen) != receipt['sequences']:
        raise ValueError('Incomplete domain input mapping')
    output.mkdir(parents=True)
    fasta = output / 'sequences.faa'
    with fasta.open('w') as handle:
        for h, sequence in sorted(unique.items()):
            handle.write(f'>S{h}\n{sequence}\n')
    links = output / 'protein_links.tsv'
    with links.open('w') as handle:
        writer = csv.DictWriter(handle, ['sequence_id'] + list(mapping[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in sorted(mapping, key=lambda r: (r['marker'], r['taxon_id'])):
            writer.writerow(dict(sequence_id='S' + row['sequence_sha256'], **row))
    result = {'source_receipt_sha256': digest(source / 'receipt.json'),
        'source_mapping_sha256': digest(source / 'protein_mapping.tsv'),
        'markers': len(files), 'taxa': len({r['taxon_id'] for r in mapping}), 'marker_proteins': len(mapping),
        'unique_sequences': len(unique), 'unique_residues': sum(len(s) for s in unique.values()),
        'deduplication': 'Exact sequence only for search; all marker, taxon and protein associations preserved.',
        'script_sha256': digest(Path(__file__)), 'artifacts': {p.name: digest(p) for p in output.iterdir()}}
    (output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    (ROOT / 'metadata/marker_domain_input_receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
