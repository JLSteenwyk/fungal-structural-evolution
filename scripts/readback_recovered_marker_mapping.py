#!/usr/bin/env python3
"""Independently reconstruct every residue link in a recovered marker snapshot."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from Bio import AlignIO, SeqIO
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mapping', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    pins = {}

    def pin(path):
        pins[str(path)] = sha(path)

    receipt_path = a.mapping / 'receipt.json'
    pin(receipt_path)
    receipt = json.loads(receipt_path.read_text())
    for name, digest in receipt['artifacts'].items():
        path = a.mapping / name
        if sha(path) != digest:
            raise ValueError('Changed mapping artifact')
        pin(path)
    links = table(a.mapping / 'marker_structure_links.tsv')
    keys = [(r['marker'], r['taxon_id']) for r in links]
    if len(keys) != len(set(keys)) or len(keys) != receipt['marker_proteins_linked']:
        raise ValueError('Link count or uniqueness differs')
    source_path = Path('results/phylogeny/markers-full-v1/protein_mapping.tsv')
    pin(source_path)
    source = {(r['marker'], r['taxon_id']): r for r in table(source_path)}
    matrix = Path('results/phylogeny/profile-matrix-50-v1')
    mr = json.loads((matrix / 'receipt.json').read_text())
    if sha(matrix / 'receipt.json') != receipt['matrix_receipt_sha256'] or sha(matrix / 'site_mapping.tsv') != mr['artifacts']['site_mapping.tsv']:
        raise ValueError('Changed matrix source')
    pin(matrix / 'receipt.json'); pin(matrix / 'site_mapping.tsv')
    sites = defaultdict(list)
    for row in table(matrix / 'site_mapping.tsv'):
        sites[row['marker']].append((int(row['alignment_column_1based']), int(row['matrix_column_1based'])))
    expected = {}
    for link in links:
        marker, taxon = link['marker'], link['taxon_id']
        original = source[marker, taxon]
        if any(link[k] != original[k] for k in ['protein_id', 'sequence_sha256']):
            raise ValueError('Marker identity differs')
        fasta = Path(f'results/phylogeny/markers-full-v1/unaligned/{marker}.faa')
        pin(fasta)
        sequences = {r.id: str(r.seq) for r in SeqIO.parse(fasta, 'fasta')}
        sequence = sequences[taxon]
        if hashlib.sha256(sequence.encode()).hexdigest() != link['sequence_sha256']:
            raise ValueError('Source sequence differs')
        prefix = Path(f'results/phylogeny/profile-alignments-full-v1/{marker}')
        pr_path, sto = prefix.with_suffix('.receipt.json'), prefix.with_suffix('.sto')
        pin(pr_path); pin(sto)
        pr = json.loads(pr_path.read_text())
        if sha(sto) != pr['stockholm_sha256']:
            raise ValueError('Alignment changed')
        aligned = {r.id: str(r.seq) for r in AlignIO.read(sto, 'stockholm')}[taxon]
        if ''.join(c for c in aligned if c not in '.-').upper() != sequence:
            raise ValueError('Alignment sequence differs')
        path = Path(link['model_path'])
        pin(path)
        if sha(path) != link['model_sha256']:
            raise ValueError('Changed coordinates')
        structure = MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False).get_structure(marker, str(path))
        residues = list(structure.get_residues())
        if ''.join(seq1(r.resname) for r in residues) != sequence:
            raise ValueError('Coordinate sequence differs')
        confidence = {r.id[1]: float(r['CA'].bfactor) for r in residues}
        if set(confidence) != set(range(1, len(sequence) + 1)):
            raise ValueError('Coordinate positions differ')
        values = []
        for profile_col, matrix_col in sites[marker]:
            column = pr['retained_stockholm_columns_1based'][profile_col - 1] - 1
            if aligned[column] in '.-':
                continue
            position = sum(c not in '.-' for c in aligned[:column + 1])
            key = marker, taxon, matrix_col
            if key in expected:
                raise ValueError('Repeated source matrix position')
            expected[key] = (link['protein_id'], link['model_id'], link['model_version'], position, confidence[position])
            values.append(confidence[position])
        if len(values) != int(link['retained_marker_residues']):
            raise ValueError('Per-link residue count differs')
        if values and (abs(sum(values) / len(values) - float(link['mean_retained_ca_plddt'])) > 1e-9 or
                       abs(sum(v >= 70 for v in values) / len(values) - float(link['fraction_retained_ca_plddt_ge70'])) > 1e-9):
            raise ValueError('Per-link confidence summary differs')
    seen = set()
    with gzip.open(a.mapping / 'matrix_to_structure_residues.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['marker'], row['taxon_id'], int(row['matrix_column_1based'])
            actual = (row['protein_id'], row['model_id'], row['model_version'], int(row['protein_residue_1based']), float(row['ca_plddt']))
            if key in seen or expected.get(key) != actual:
                raise ValueError('Unexpected, duplicate, or incorrect residue link')
            seen.add(key)
    if seen != set(expected) or len(seen) != receipt['matrix_residue_links']:
        raise ValueError('Incomplete residue mapping')
    for path, digest in pins.items():
        if sha(path) != digest:
            raise ValueError('Input changed during audit')
    result = dict(status='passed_complete_recovered_marker_residue_readback', links=len(links),
                  residue_links=len(seen), script_sha256=sha(__file__), source_sha256=pins,
                  scope='Every emitted link and omitted alignment gap checked by prefix-count projection; coordinate-derived sequences, CA confidence, and per-link summaries verified. Does not merge or fit phylogenetic datasets.')
    with a.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}))


if __name__ == '__main__':
    main()
