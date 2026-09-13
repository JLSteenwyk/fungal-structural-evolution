#!/usr/bin/env python3
"""Validate native 3Di descriptors and record confidence of all six feature residues."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT, sha
from retrieve_marker_pae import validate_pae


def normalize(v):
    lengths = np.linalg.norm(v, axis=-1, keepdims=True)
    if np.any(lengths == 0) or not np.isfinite(lengths).all():
        raise ValueError('Degenerate or nonfinite coordinate vectors')
    return v / lengths


def coordinate_features(ca, n, c, cb):
    """Reproduce pinned Foldseek partner selection/features, not neural encoding."""
    ca, n, c, cb = (np.asarray(x, float).copy() for x in [ca, n, c, cb])
    length = len(ca)
    if length < 4 or any(x.shape != (length, 3) for x in [ca, n, c, cb]):
        raise ValueError('Complete single-chain backbone arrays required')
    if not all(np.isfinite(x).all() for x in [ca, n, c]):
        raise ValueError('Incomplete backbone requires separate validity review')
    missing = np.isnan(cb[:, 0])
    if missing.any():
        v1, v2 = normalize(c[missing] - ca[missing]), normalize(n[missing] - ca[missing])
        b1 = v2 + v1 / 3
        u1, u2 = normalize(b1), normalize(np.cross(v1, b1))
        v4 = -v1 / 3 + (-u1 / 2 - u2 * np.sqrt(3) / 2) * np.sqrt(8) / 3
        cb[missing] = ca[missing] + v4 * 1.5336
    v = cb - ca
    k = normalize(np.cross(v, n - ca))
    alpha = 270 / 180 * 3.14159265359
    rotated = v * np.cos(alpha) + np.cross(k, v) * np.sin(alpha) + k * np.sum(k * v, axis=1)[:, None] * (1 - np.cos(alpha))
    # The pinned second rotation has beta=0 and the center-distance factor is 2.
    centers = ca + rotated * 2
    distances = np.linalg.norm(centers[:, None] - centers[None, :], axis=2)
    distances[:, [0, length - 1]] = np.inf
    np.fill_diagonal(distances, np.inf)
    partner = np.argmin(distances, axis=1)
    valid = np.ones(length, dtype=bool)
    valid[[0, length - 1]] = False
    i = np.arange(1, length - 1)
    j = partner[i]
    u1, u2 = normalize(ca[i] - ca[i - 1]), normalize(ca[i + 1] - ca[i])
    u3, u4 = normalize(ca[j] - ca[j - 1]), normalize(ca[j + 1] - ca[j])
    u5 = normalize(ca[j] - ca[i])
    features = np.zeros((length, 10))
    for column, (a, b) in enumerate([(u1, u2), (u3, u4), (u1, u5), (u3, u5), (u1, u4), (u2, u3), (u1, u3)]):
        features[i, column] = np.sum(a * b, axis=1)
    features[i, 7] = np.linalg.norm(ca[i] - ca[j], axis=1)
    features[i, 8] = np.clip(j - i, -4, 4)
    features[i, 9] = np.sign(j - i) * np.log(np.abs(j - i) + 1)
    partner[~valid] = -1
    return partner, valid, features


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--pae', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output = args.output.resolve()
    checked_receipt(args.snapshot)
    pr = checked_receipt(args.pae)
    config = json.loads((args.native / 'config.json').read_text())
    if config['mapping_receipt_sha256'] != sha(args.snapshot / 'receipt.json') or pr['mapping_receipt_sha256'] != config['mapping_receipt_sha256']:
        raise ValueError('Native extraction, coordinates and PAE have different provenance')
    for name, checksum in config['source_sha256'].items():
        if sha(args.native / 'source' / name) != checksum:
            raise ValueError('Changed pinned Foldseek source')
    models = json.loads((args.snapshot / 'model_provenance.json').read_text())
    if models != json.loads((args.native / 'model_provenance.json').read_text()):
        raise ValueError('Native extraction model inventory differs')
    by_name = {Path(r['path']).stem: r for r in models}
    aa_rows = list(SeqIO.parse(args.native / 'amino_acids.faa', 'fasta'))
    ss_rows = list(SeqIO.parse(args.native / 'states_3di.faa', 'fasta'))
    aas, states = {r.id: str(r.seq) for r in aa_rows}, {r.id: str(r.seq) for r in ss_rows}
    if len(aas) != len(aa_rows) or len(states) != len(ss_rows) or set(aas) != set(by_name) or set(states) != set(by_name):
        raise ValueError('Native FASTA identity/count mismatch')
    descriptors = {}
    with (args.native / 'descriptors.tsv').open() as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 4:
                raise ValueError('Unexpected native descriptor schema')
            name = Path(fields[0].split()[0]).stem
            if name in descriptors or name not in by_name or fields[1] != aas[name] or fields[2] != states[name]:
                raise ValueError('Native descriptor and DB exports disagree')
            data = np.array([float(x) for x in fields[3].split(',')]).reshape(-1, 10)
            if data.shape != (len(aas[name]), 10) or not np.isfinite(data).all():
                raise ValueError('Invalid native feature array')
            descriptors[name] = data
    if set(descriptors) != set(by_name):
        raise ValueError('Incomplete descriptor model coverage')
    pae_rows = {(r['model_id'], str(r['version'])): r for r in json.loads((args.pae / 'pae_manifest.json').read_text()) if r['status'] == 'verified'}
    if args.output.exists():
        raise FileExistsError('Use a new immutable audited encoding directory')
    args.output.mkdir(parents=True)
    summaries = []
    for name, model in sorted(by_name.items()):
        sequence, structural = aas[name], states[name]
        if hashlib.sha256(sequence.encode()).hexdigest() != model['sequence_sha256'] or len(structural) != len(sequence) or not set(structural) <= set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('Native amino-acid sequence identity or 3Di alphabet differs')
        path = ROOT / model['path']
        if sha(path) != model['sha256']:
            raise ValueError('Changed coordinates')
        cif = MMCIF2Dict(str(path))
        if [''.join(s.split()) for s in cif['_entity_poly.pdbx_seq_one_letter_code_can']] != [sequence]:
            raise ValueError('Coordinate and exported complete sequences disagree')
        length = len(sequence)
        atoms = {a: np.full((length, 3), np.nan) for a in ['CA', 'N', 'C', 'CB']}
        plddt = np.full(length, np.nan)
        seen = set()
        for atom, pos, x, y, z, confidence in zip(cif['_atom_site.label_atom_id'], cif['_atom_site.label_seq_id'],
                cif['_atom_site.Cartn_x'], cif['_atom_site.Cartn_y'], cif['_atom_site.Cartn_z'], cif['_atom_site.B_iso_or_equiv']):
            if atom not in atoms:
                continue
            i = int(pos) - 1
            if not 0 <= i < length or (atom, i) in seen:
                raise ValueError('Ambiguous atom/residue identity')
            seen.add((atom, i)); atoms[atom][i] = [float(x), float(y), float(z)]
            if atom == 'CA':
                plddt[i] = float(confidence)
        if not np.isfinite(plddt).all() or np.any(plddt < 0) or np.any(plddt > 100):
            raise ValueError('Invalid CA confidence')
        partner, valid, features = coordinate_features(atoms['CA'], atoms['N'], atoms['C'], atoms['CB'])
        # Native SSTR exports four significant digits; partner recovery by log inversion alone can be ambiguous.
        if not np.array_equal(valid, descriptors[name][:, 9] != 0) or not np.allclose(features, descriptors[name], rtol=5.2e-4, atol=1e-10):
            error = np.abs(features - descriptors[name])
            raise ValueError(f'Native feature/partner reconstruction failed for {name}: max discrepancy {error.max()}')
        pae_row = pae_rows[(model['model_id'], str(model['version']))]
        pae_path = ROOT / pae_row['path']
        if pae_row['sequence_sha256'] != model['sequence_sha256'] or sha(pae_path) != pae_row['gzip_sha256']:
            raise ValueError('PAE identity differs')
        pae = validate_pae(gzip.decompress(pae_path.read_bytes()), length)
        i = np.flatnonzero(valid); j = partner[i]
        context = np.stack([i - 1, i, i + 1, j - 1, j, j + 1], axis=1)
        minimum, maximum = np.full(length, np.nan), np.full(length, np.nan)
        minimum[i] = plddt[context].min(axis=1)
        maximum[i] = pae[context[:, :, None], context[:, None, :]].max(axis=(1, 2))
        output = args.output / (name + '.npz')
        np.savez_compressed(output, sequence=np.array(sequence), states=np.array(structural), valid=valid,
            partner_residue_1based=np.where(valid, partner + 1, 0), ca_plddt=plddt,
            feature_min_plddt=minimum, feature_max_pae=maximum)
        summaries.append({'model_name': name, 'model_id': model['model_id'], 'version': model['version'],
            'sequence_sha256': model['sequence_sha256'], 'length': length, 'valid_states': int(valid.sum()),
            'invalid_states': int((~valid).sum()), 'valid_focal_plddt70': int((valid & (plddt >= 70)).sum()),
            'valid_feature_plddt70': int((valid & (minimum >= 70)).sum()),
            'valid_feature_plddt70_pae10': int((valid & (minimum >= 70) & (maximum <= 10)).sum()),
            'encoding_path': str(output.relative_to(ROOT)), 'encoding_sha256': sha(output)})
    path = args.output / 'model_summary.tsv'
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, list(summaries[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(summaries)
    receipt = {'status': 'complete_native_3di_feature_audit', 'models': len(summaries),
        'totals': {k: sum(r[k] for r in summaries) for k in ['length', 'valid_states', 'invalid_states', 'valid_focal_plddt70', 'valid_feature_plddt70', 'valid_feature_plddt70_pae10']},
        'native_config_sha256': sha(args.native / 'config.json'), 'native_path': str(args.native),
        'native_artifacts': {n: sha(args.native / n) for n in ['amino_acids.faa', 'states_3di.faa', 'descriptors.tsv', 'model_provenance.json']},
        'mapping_receipt_sha256': sha(args.snapshot / 'receipt.json'), 'pae_receipt_sha256': sha(args.pae / 'receipt.json'),
        'script_sha256': sha(Path(__file__)), 'numpy_version': np.__version__,
        'feature_verification_tolerance': {'rtol': 5.2e-4, 'atol': 1e-10, 'reason': 'Native descriptor output has four significant digits'},
        'interpretation': 'Coordinate-derived 20-state alphabet with explicit validity mask. Invalid terminal states use the ordinary coil symbol and must not be counted as observations. Confidence contexts cover i-1,i,i+1,j-1,j,j+1; PAE is the maximum over all directional pairs in that context. Alphabet states are not amino acids and distances are not physical displacement.',
        'artifacts': {'model_summary.tsv': sha(path)}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
