#!/usr/bin/env python3
"""Assess matched residue-distance changes under directional PAE confidence filters."""
import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from compare_marker_structures import AA, ROOT, geometry, sha
from retrieve_marker_pae import validate_pae


def confidence_mask(pae_a, pae_b, positions_a, positions_b, cutoff):
    """Require both directional errors in both models, using 1-based residue IDs."""
    a, b = np.asarray(positions_a, dtype=int) - 1, np.asarray(positions_b, dtype=int) - 1
    if min(a.min(), b.min()) < 0 or a.max() >= len(pae_a) or b.max() >= len(pae_b):
        raise ValueError('Residue outside PAE matrix')
    x, y = pae_a[np.ix_(a, a)], pae_b[np.ix_(b, b)]
    return np.maximum.reduce([x, x.T, y, y.T]) <= cutoff


def checked_receipt(directory):
    receipt = json.loads((directory / 'receipt.json').read_text())
    for filename, digest in receipt['artifacts'].items():
        if sha(directory / filename) != digest:
            raise ValueError(f'Changed artifact: {directory / filename}')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--comparisons', type=Path, required=True)
    parser.add_argument('--pae', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    mapping_receipt = checked_receipt(args.snapshot)
    comparisons_receipt = checked_receipt(args.comparisons)
    pae_receipt = checked_receipt(args.pae)
    for receipt in [comparisons_receipt, pae_receipt]:
        if receipt['mapping_receipt_sha256'] != sha(args.snapshot / 'receipt.json'):
            raise ValueError('Inputs refer to different mapping snapshots')
    if args.output.exists():
        raise FileExistsError('Use a new immutable output directory')
    matrix_dir = ROOT / 'results/phylogeny/profile-matrix-50-v1'
    if sha(matrix_dir / 'receipt.json') != mapping_receipt['matrix_receipt_sha256']:
        raise ValueError('Changed species matrix receipt')
    checked_receipt(matrix_dir)
    matrix = {r.id: str(r.seq) for r in SeqIO.parse(matrix_dir / 'matrix.faa', 'fasta')}
    with (args.snapshot / 'marker_structure_links.tsv').open() as handle:
        links = {(r['marker'], r['taxon_id']): r for r in csv.DictReader(handle, delimiter='\t')}
    mapped = defaultdict(dict)
    with gzip.open(args.snapshot / 'matrix_to_structure_residues.tsv.gz', 'rt') as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            mapped[(r['marker'], r['taxon_id'])][int(r['matrix_column_1based'])] = (
                int(r['protein_residue_1based']), float(r['ca_plddt']))
    paes = {}
    for r in json.loads((args.pae / 'pae_manifest.json').read_text()):
        if r['status'] != 'verified':
            continue
        path = ROOT / r['path']
        if sha(path) != r['gzip_sha256']:
            raise ValueError('Changed PAE matrix')
        paes[(r['model_id'], str(r['version']))] = (validate_pae(gzip.decompress(path.read_bytes()), r['length']), r)
    coordinates = {}
    for r in links.values():
        if r['model_path'] in coordinates:
            continue
        path = ROOT / r['model_path']
        if sha(path) != r['model_sha256']:
            raise ValueError('Changed model coordinates')
        cif = MMCIF2Dict(str(path))
        coordinates[r['model_path']] = {int(p): np.array([float(x), float(y), float(z)])
            for atom, p, x, y, z in zip(cif['_atom_site.label_atom_id'], cif['_atom_site.label_seq_id'],
                cif['_atom_site.Cartn_x'], cif['_atom_site.Cartn_y'], cif['_atom_site.Cartn_z']) if atom == 'CA'}
    with (args.comparisons / 'pairwise_metrics.tsv').open() as handle:
        pairs = [r for r in csv.DictReader(handle, delimiter='\t') if int(r['plddt_cutoff']) == 70]
    output, exclusions = [], []
    for baseline in pairs:
        marker = baseline['marker']
        ka, kb = (marker, baseline['taxon_a']), (marker, baseline['taxon_b'])
        a, b = links[ka], links[kb]
        pa, pb = paes.get((a['model_id'], a['model_version'])), paes.get((b['model_id'], b['model_version']))
        base = {k: baseline[k] for k in ['marker', 'taxon_a', 'taxon_b', 'same_model_coordinates',
                 'compared_residues', 'uncorrected_sequence_difference', 'ca_superposition_rmsd_angstrom']}
        if pa is None or pb is None:
            exclusions.append(dict(base, reason='Missing verified PAE for one or both models'))
            continue
        for link, pae in [(a, pa), (b, pb)]:
            if link['sequence_sha256'] != pae[1]['sequence_sha256']:
                raise ValueError('PAE and coordinate sequence provenance differ')
        ma, mb = mapped[ka], mapped[kb]
        columns = [c for c in sorted(set(ma) & set(mb)) if ma[c][1] >= 70 and mb[c][1] >= 70
                   and matrix[ka[1]][c - 1] in AA and matrix[kb[1]][c - 1] in AA]
        pos_a, pos_b = np.array([ma[c][0] for c in columns]), np.array([mb[c][0] for c in columns])
        x = np.array([coordinates[a['model_path']][p] for p in pos_a])
        y = np.array([coordinates[b['model_path']][p] for p in pos_b])
        if len(columns) != int(baseline['compared_residues']):
            raise ValueError('Qualified residue set differs from direct comparison')
        g = geometry(x, y, pos_a, pos_b)
        if not np.isclose(g['ca_superposition_rmsd_angstrom'], float(baseline['ca_superposition_rmsd_angstrom']), atol=1e-8):
            raise ValueError('Baseline geometry not reproduced')
        dx = np.linalg.norm(x[:, None] - x[None, :], axis=2)
        dy = np.linalg.norm(y[:, None] - y[None, :], axis=2)
        eligible = (np.triu(np.ones(dx.shape, dtype=bool), 1)
                    & (np.abs(pos_a[:, None] - pos_a[None, :]) >= 3)
                    & (np.abs(pos_b[:, None] - pos_b[None, :]) >= 3))
        local = eligible & ((dx <= 15) | (dy <= 15))
        if int(local.sum()) != int(baseline['local_distance_pairs']):
            raise ValueError('Baseline local residue pairs not reproduced')
        delta = np.abs(dx - dy)
        for cutoff in [5, 10, 15]:
            confident = confidence_mask(pa[0], pb[0], pos_a, pos_b, cutoff)
            for scope, mask in [('all_nonadjacent', eligible), ('local_15A_union', local)]:
                good, uncertain = mask & confident, mask & ~confident
                output.append(dict(base, plddt_cutoff=70, pae_cutoff_angstrom=cutoff,
                    residue_pair_scope=scope, eligible_residue_pairs=int(mask.sum()),
                    confident_residue_pairs=int(good.sum()), uncertain_residue_pairs=int(uncertain.sum()),
                    confident_fraction=float(good.sum() / mask.sum()) if mask.any() else '',
                    unfiltered_mean_absolute_distance_change_angstrom=float(delta[mask].mean()) if mask.any() else '',
                    confident_mean_absolute_distance_change_angstrom=float(delta[good].mean()) if good.any() else '',
                    uncertain_mean_absolute_distance_change_angstrom=float(delta[uncertain].mean()) if uncertain.any() else ''))
    args.output.mkdir(parents=True)
    for filename, rows in [('pae_sensitivity.tsv', output), ('excluded_pairs.tsv', exclusions)]:
        with (args.output / filename).open('w') as handle:
            if rows:
                writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows(rows)
    receipt = {'mapping_receipt_sha256': sha(args.snapshot / 'receipt.json'),
        'comparisons_receipt_sha256': sha(args.comparisons / 'receipt.json'),
        'pae_receipt_sha256': sha(args.pae / 'receipt.json'), 'script_sha256': sha(Path(__file__)),
        'requested_taxon_pairs': len(pairs), 'assessed_taxon_pairs': len(pairs) - len(exclusions),
        'excluded_taxon_pairs': len(exclusions), 'sensitivity_rows': len(output),
        'interpretation': 'Descriptive confidence sensitivity, not significance tests, domain assignments or branch rates. PAE is predicted alignment uncertainty, not a calibrated error bar on residue distances. Filters change residue-pair composition.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
