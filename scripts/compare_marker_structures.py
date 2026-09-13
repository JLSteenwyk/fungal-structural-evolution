#!/usr/bin/env python3
"""Compute descriptive, alignment-matched CA geometry and sequence differences."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

ROOT = Path(__file__).resolve().parents[1]
AA = set('ACDEFGHIKLMNPQRSTVWY')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geometry(x, y, positions_x, positions_y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 2 or x.shape[1] != 3 or len(x) < 3:
        raise ValueError('At least three matched three-dimensional coordinates required')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Nonfinite coordinates')
    xc, yc = x - x.mean(axis=0), y - y.mean(axis=0)
    u, _, vt = np.linalg.svd(xc.T @ yc)
    correction = np.diag([1., 1., 1. if np.linalg.det(u @ vt) >= 0 else -1.])
    rotation = u @ correction @ vt
    rmsd = np.sqrt(np.mean(np.sum((xc @ rotation - yc) ** 2, axis=1)))
    dx = np.linalg.norm(x[:, None] - x[None, :], axis=2)
    dy = np.linalg.norm(y[:, None] - y[None, :], axis=2)
    px, py = np.asarray(positions_x), np.asarray(positions_y)
    neighbors = (np.triu(np.ones(dx.shape, dtype=bool), 1)
                 & ((dx <= 15) | (dy <= 15))
                 & (np.abs(px[:, None] - px[None, :]) >= 3)
                 & (np.abs(py[:, None] - py[None, :]) >= 3))
    differences = dx[neighbors] - dy[neighbors]
    return {'ca_superposition_rmsd_angstrom': float(rmsd),
            'local_distance_pairs': int(len(differences)),
            'local_distance_mean_absolute_change_angstrom': float(np.mean(np.abs(differences))) if len(differences) else '',
            'local_distance_rms_change_angstrom': float(np.sqrt(np.mean(differences ** 2))) if len(differences) else ''}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads((args.snapshot / 'receipt.json').read_text())
    for filename, checksum in receipt['artifacts'].items():
        if sha(args.snapshot / filename) != checksum:
            raise ValueError('Changed structural mapping artifact')
    if args.output.exists():
        raise FileExistsError('Use an immutable new comparison snapshot')
    with (args.snapshot / 'marker_structure_links.tsv').open() as handle:
        links = list(csv.DictReader(handle, delimiter='\t'))
    groups = defaultdict(list)
    for row in links:
        groups[row['marker']].append(row)
    mapped = defaultdict(dict)
    with gzip.open(args.snapshot / 'matrix_to_structure_residues.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            mapped[(row['marker'], row['taxon_id'])][int(row['matrix_column_1based'])] = (
                int(row['protein_residue_1based']), float(row['ca_plddt']))
    matrix_dir = ROOT / 'results/phylogeny/profile-matrix-50-v1'
    if sha(matrix_dir / 'receipt.json') != receipt['matrix_receipt_sha256']:
        raise ValueError('Changed matrix receipt')
    matrix_receipt = json.loads((matrix_dir / 'receipt.json').read_text())
    if sha(matrix_dir / 'matrix.faa') != matrix_receipt['artifacts']['matrix.faa']:
        raise ValueError('Changed species matrix')
    with (matrix_dir / 'matrix.faa').open() as handle:
        matrix = {r.id: str(r.seq) for r in SeqIO.parse(handle, 'fasta')}
    models = {}
    for row in links:
        path = ROOT / row['model_path']
        if row['model_path'] in models:
            continue
        if sha(path) != row['model_sha256']:
            raise ValueError('Changed model coordinates')
        cif = MMCIF2Dict(str(path))
        models[row['model_path']] = {int(pos): np.array([float(x), float(y), float(z)])
            for atom, pos, x, y, z in zip(cif['_atom_site.label_atom_id'], cif['_atom_site.label_seq_id'],
                cif['_atom_site.Cartn_x'], cif['_atom_site.Cartn_y'], cif['_atom_site.Cartn_z']) if atom == 'CA'}
    args.output.mkdir(parents=True)
    results, exclusions = [], []
    for marker, rows in sorted(groups.items()):
        for a, b in itertools.combinations(sorted(rows, key=lambda r: r['taxon_id']), 2):
            ma, mb = mapped[(marker, a['taxon_id'])], mapped[(marker, b['taxon_id'])]
            common = sorted(set(ma) & set(mb))
            for cutoff in [50, 70, 90]:
                columns = [c for c in common if ma[c][1] >= cutoff and mb[c][1] >= cutoff
                           and matrix[a['taxon_id']][c - 1] in AA and matrix[b['taxon_id']][c - 1] in AA]
                base = {'marker': marker, 'taxon_a': a['taxon_id'], 'taxon_b': b['taxon_id'],
                        'plddt_cutoff': cutoff, 'shared_profile_positions': len(common), 'compared_residues': len(columns),
                        'fraction_of_shared_positions_compared': len(columns) / len(common) if common else 0,
                        'same_model_coordinates': a['model_sha256'] == b['model_sha256']}
                if len(columns) < 50 or len(columns) < .5 * len(common):
                    exclusions.append(dict(base, reason='Fewer than 50 qualified residues or less than half shared positions'))
                    continue
                pa, pb = [ma[c][0] for c in columns], [mb[c][0] for c in columns]
                x = [models[a['model_path']][p] for p in pa]
                y = [models[b['model_path']][p] for p in pb]
                difference = sum(matrix[a['taxon_id']][c - 1] != matrix[b['taxon_id']][c - 1] for c in columns) / len(columns)
                results.append(dict(base, uncorrected_sequence_difference=difference, **geometry(x, y, pa, pb)))
    for filename, rows in [('pairwise_metrics.tsv', results), ('pairwise_exclusions.tsv', exclusions)]:
        with (args.output / filename).open('w') as handle:
            if rows:
                writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows(rows)
    result = {'mapping_receipt_sha256': sha(args.snapshot / 'receipt.json'), 'script_sha256': sha(Path(__file__)), 'numpy_version': np.__version__,
              'comparison_rows': len(results), 'excluded_rows': len(exclusions),
              'distinct_marker_taxon_pairs': len({(r['marker'], r['taxon_a'], r['taxon_b']) for r in results}),
              'markers_compared': len({r['marker'] for r in results}),
              'interpretation': 'Descriptive matched-site metrics, not branch lengths, rates, independent observations or evidence of selection. Domain motion and prediction artifacts remain possible.',
              'local_metric_definition': 'CA distance changes for unordered residue pairs within 15 Angstrom in either model, separated by at least three sequence positions in both. This is not standard lDDT.',
              'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
