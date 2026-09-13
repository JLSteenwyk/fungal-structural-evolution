#!/usr/bin/env python3
"""Read back native sequence/state/descriptor exports before PAE-dependent auditing."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha


def fasta(path):
    rows = list(SeqIO.parse(path, 'fasta'))
    result = {r.id: str(r.seq) for r in rows}
    if len(result) != len(rows):
        raise ValueError('Duplicate FASTA model identity')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new native-export audit receipt')
    checked_receipt(args.snapshot)
    config = json.loads((args.native / 'config.json').read_text())
    if config['mapping_receipt_sha256'] != sha(args.snapshot / 'receipt.json'):
        raise ValueError('Mapping provenance differs')
    for name, checksum in config['source_sha256'].items():
        if sha(args.native / 'source' / name) != checksum:
            raise ValueError('Changed pinned native source')
    models = json.loads((args.snapshot / 'model_provenance.json').read_text())
    if models != json.loads((args.native / 'model_provenance.json').read_text()):
        raise ValueError('Native model inventory differs')
    by_name = {Path(r['path']).stem: r for r in models}
    aas = fasta(args.native / 'amino_acids.faa')
    states = fasta(args.native / 'states_3di.faa')
    if len(by_name) != len(models) or set(aas) != set(by_name) or set(states) != set(by_name):
        raise ValueError('Incomplete or ambiguous native model coverage')
    for name, model in by_name.items():
        if hashlib.sha256(aas[name].encode()).hexdigest() != model['sequence_sha256'] or len(aas[name]) != model['length']:
            raise ValueError('Native complete protein sequence differs')
        if len(states[name]) != model['length'] or not set(states[name]) <= set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('Invalid state alphabet or sequence length')
    seen = set()
    residues = zero_offset = 0
    with (args.native / 'descriptors.tsv').open() as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 4:
                raise ValueError('Unexpected descriptor schema')
            name = Path(fields[0].split()[0]).stem
            if name in seen or name not in by_name or fields[1] != aas[name] or fields[2] != states[name]:
                raise ValueError('Descriptor exports disagree or repeat a model')
            data = np.array([float(x) for x in fields[3].split(',')]).reshape(-1, 10)
            if data.shape != (len(aas[name]), 10) or not np.isfinite(data).all():
                raise ValueError('Incomplete/nonfinite native descriptor')
            residues += len(data)
            zero_offset += int((data[:, 9] == 0).sum())
            seen.add(name)
    if seen != set(by_name):
        raise ValueError('Incomplete native descriptor model set')
    result = {'status': 'complete_native_export_identity_audit', 'models': len(models),
        'residues': residues, 'descriptor_values': residues * 10,
        'zero_log_sequence_offset_rows': zero_offset,
        'native_path': str(args.native), 'mapping_receipt_sha256': sha(args.snapshot / 'receipt.json'),
        'native_config_sha256': sha(args.native / 'config.json'), 'script_sha256': sha(Path(__file__)),
        'interpretation': 'Exact full-sequence identity, state length/alphabet, export agreement and finite ten-feature shape checked for every model. Zero-offset descriptor rows are counted without accepting their state as valid. Coordinate-derived feature/partner reconstruction and six-residue pLDDT/PAE qualification are still required; this receipt does not authorize evolutionary use of unmasked states.',
        'artifacts': {name: sha(args.native / name) for name in ['config.json', 'model_provenance.json', 'amino_acids.faa', 'states_3di.faa', 'descriptors.tsv']}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
