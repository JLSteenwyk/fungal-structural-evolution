#!/usr/bin/env python3
"""Independently read back completed prediction artifacts and summarize coverage/resources."""
import argparse
import csv
import io
import json
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from prepare_pfam import digest


def normalize_link(row):
    row = dict(row)
    sid = row.get('sequence_id')
    digest_value = row.get('sequence_sha256')
    if sid and digest_value and sid != 'S' + digest_value:
        raise ValueError('Conflicting link sequence identity')
    if not sid and digest_value:
        sid = 'S' + digest_value
    if not sid or len(sid) != 65 or sid[0] != 'S' or any(c not in '0123456789abcdef' for c in sid[1:]):
        raise ValueError('Missing or malformed link sequence identity')
    row['sequence_id'] = sid
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predictions', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--links', type=Path, help='Explicit marker/reference links; default inputs/all_marker_links.tsv')
    args = parser.parse_args()
    chunk = json.loads((args.predictions / 'last_chunk.json').read_text())
    config_path = args.predictions / 'config.json'
    config = json.loads(config_path.read_text())
    if chunk['config_sha256'] != digest(config_path) or config['input_receipt_sha256'] != digest(args.inputs / 'receipt.json'):
        raise ValueError('Changed source configuration')
    ir = json.loads((args.inputs / 'receipt.json').read_text())
    for name, sha in ir['artifacts'].items():
        if digest(args.inputs / name) != sha:
            raise ValueError('Changed prediction inputs')
    sequences = {r.id: str(r.seq) for r in SeqIO.parse(args.inputs / 'candidates.faa', 'fasta')}
    link_path = args.links or args.inputs / 'all_marker_links.tsv'
    links = [normalize_link(r) for r in csv.DictReader(link_path.open(), delimiter='\t')]
    records = []
    pdb_parser = PDBParser(QUIET=True)
    for path in sorted(args.predictions.glob('S*.json')):
        row = json.loads(path.read_text())
        if row['status'] != 'verified_prediction':
            continue
        sid = row['sequence_id']
        if row['config_sha256'] != digest(config_path) or sid not in sequences:
            raise ValueError('Prediction identity/configuration mismatch')
        for name, sha in row['artifacts'].items():
            if digest(args.predictions / name) != sha:
                raise ValueError('Changed prediction artifact')
        sequence = sequences[sid]
        model = pdb_parser.get_structure(sid, io.StringIO((args.predictions / (sid + '.pdb')).read_text()))
        residues = list(model.get_residues())
        if (len(list(model.get_models())) != 1 or len(list(model.get_chains())) != 1
                or ''.join(seq1(r.resname) for r in residues) != sequence
                or [r.id[1] for r in residues] != list(range(1, len(sequence) + 1))
                or any('CA' not in r for r in residues)):
            raise ValueError('PDB readback failed sequence/numbering/chain check')
        if not all(np.isfinite(atom.coord).all() for atom in model.get_atoms()):
            raise ValueError('Nonfinite atom coordinate')
        with np.load(args.predictions / (sid + '.npz'), allow_pickle=False) as data:
            plddt, pae = data['ca_plddt'], data['pae']
            if (str(data['sequence']) != sequence or plddt.shape != (len(sequence),)
                    or pae.shape != (len(sequence), len(sequence))
                    or not np.isfinite(plddt).all() or not np.isfinite(pae).all()
                    or np.any(plddt < 0) or np.any(plddt > 100) or np.any(pae < 0)
                    or np.any(pae > float(data['max_pae']) + 1e-4)
                    or not np.allclose([r['CA'].bfactor for r in residues], plddt, atol=.0051, rtol=0)
                    or abs(float(plddt.mean()) - row['mean_ca_plddt']) > 1e-5):
                raise ValueError('Confidence/PAE readback validation failed')
        records.append({k: row[k] for k in ['sequence_id', 'length', 'mean_ca_plddt',
                         'fraction_ca_plddt_below50', 'inference_seconds', 'peak_gpu_allocated_bytes']}
                       | {'prediction_receipt_sha256': digest(path)})
    ids = {r['sequence_id'] for r in records}
    if len(ids) != len(records) or len(records) != chunk['cached_predictions'] + chunk['new_predictions']:
        raise ValueError('Completed chunk and actual predictions disagree')
    matched_links = [r for r in links if r['sequence_id'] in ids]
    if args.output.exists():
        raise FileExistsError('Use an immutable new audit snapshot')
    args.output.mkdir(parents=True)
    for name, rows in [('predictions.tsv', records), ('taxon_links.tsv', matched_links)]:
        with (args.output / name).open('w') as handle:
            if rows:
                writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows(rows)
    summary = {'status': 'complete_artifact_readback', 'predictions': len(records),
        'taxon_marker_links': len(matched_links), 'taxa': len({r['taxon_id'] for r in matched_links}),
        'markers': len({r['marker'] for r in matched_links}),
        'links_path': str(link_path), 'links_sha256': digest(link_path),
        'config_sha256': digest(config_path), 'chunk_receipt_sha256': digest(args.predictions / 'last_chunk.json'),
        'script_sha256': digest(Path(__file__)), 'remaining_eligible': chunk['remaining_eligible'],
        'length_or_alphabet_deferred': chunk['length_or_alphabet_deferred'],
        'interpretation': 'Short-protein production chunk; verified serialization and confidence provenance, not experimental accuracy or completed structural atlas.',
        'artifacts': {p.name: digest(p) for p in args.output.iterdir()}}
    for field in ['length', 'inference_seconds', 'peak_gpu_allocated_bytes', 'mean_ca_plddt']:
        values = [r[field] for r in records]
        summary[field] = {'min': min(values), 'median': float(np.median(values)), 'max': max(values)} if values else None
    (args.output / 'receipt.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
