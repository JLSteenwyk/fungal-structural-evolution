#!/usr/bin/env python3
"""Validate recovered marker candidates using coordinate-derived sequences."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--screen', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    receipt_path = args.screen / 'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    candidates_path = args.screen / 'candidate_models.jsonl'
    digest = sha(candidates_path)
    if digest != receipt['artifacts']['candidate_models.jsonl']:
        raise ValueError('Candidate inventory changed')
    candidates = [json.loads(line) for line in candidates_path.read_text().splitlines()]
    if len(candidates) != receipt['candidate_models']:
        raise ValueError('Candidate count differs')
    results = []
    for entry in candidates:
        model = entry['model']
        path = Path(model['path'])
        if sha(path) != model['sha256']:
            raise ValueError('Coordinate bytes changed')
        structure = MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False).get_structure(model['model_id'], str(path))
        models = list(structure)
        if len(models) != 1 or len(list(models[0])) != 1:
            raise ValueError('Expected one model and one chain')
        residues = list(next(iter(models[0])))
        if [r.id for r in residues] != [(' ', i, ' ') for i in range(1, model['length'] + 1)]:
            raise ValueError('Incomplete or unexpected residue numbering')
        sequence = ''.join(seq1(r.resname) for r in residues)
        if hashlib.sha256(sequence.encode()).hexdigest() != model['sequence_sha256']:
            raise ValueError('Coordinate-derived sequence differs')
        confidence = []
        for residue in residues:
            if residue.is_disordered() or 'CA' not in residue:
                raise ValueError('Missing or ambiguous residue')
            for atom in residue:
                if atom.is_disordered() or not all(math.isfinite(float(x)) for x in atom.coord):
                    raise ValueError('Ambiguous or nonfinite atom')
            value = float(residue['CA'].bfactor)
            if not math.isfinite(value) or not 0 <= value <= 100:
                raise ValueError('Invalid confidence')
            confidence.append(value)
        mean = sum(confidence) / len(confidence)
        below50 = sum(v < 50 for v in confidence) / len(confidence)
        if abs(mean - model['mean_ca_plddt']) > 1e-6 or abs(below50 - model['fraction_ca_plddt_below50']) > 1e-9:
            raise ValueError('Stored confidence summary differs')
        if sha(path) != model['sha256']:
            raise ValueError('Coordinates changed during readback')
        results.append(dict(model_id=model['model_id'], version=model['version'],
                            path=str(path), coordinate_sha256=model['sha256'],
                            sequence_sha256=model['sequence_sha256'], length=len(sequence),
                            mean_ca_plddt=mean, ca_plddt_below50=sum(v < 50 for v in confidence),
                            ca_plddt_at_least70=sum(v >= 70 for v in confidence),
                            ca_plddt_at_least90=sum(v >= 90 for v in confidence)))
    output = dict(status='passed_coordinate_sequence_and_confidence_readback',
                  screen_receipt_sha256=sha(receipt_path), candidate_inventory_sha256=digest,
                  script_sha256=sha(__file__), models=len(results),
                  residues=sum(r['length'] for r in results), results=results,
                  scope='Independent MMCIFParser coordinate traversal; exact full-sequence hash, complete CA coverage, finite coordinates, and residue confidence summaries. Not PAE qualification, structural-alphabet encoding, or integration into frozen analyses.')
    with args.output.open('x') as handle:
        handle.write(json.dumps(output, indent=2) + '\n')
    print(json.dumps({k: v for k, v in output.items() if k != 'results'}, indent=2))


if __name__ == '__main__':
    main()
