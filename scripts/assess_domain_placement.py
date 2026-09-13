#!/usr/bin/env python3
"""Assess confidence and residue-distance differences between matched domain instances."""
import argparse
import csv
import gzip
import itertools
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from compare_marker_structures import AA, ROOT, sha
from assess_pae_sensitivity import checked_receipt, confidence_mask
from retrieve_marker_pae import validate_pae


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparisons', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt = checked_receipt(args.comparisons)
    if args.output.exists():
        raise FileExistsError('Use an immutable output snapshot')
    sources = {k: Path(v['path']) for k, v in receipt['sources'].items()}
    for key in ['snapshot', 'pae']:
        checked_receipt(sources[key])
        if sha(sources[key] / 'receipt.json') != receipt['sources'][key]['receipt_sha256']:
            raise ValueError('Changed source receipt')
    with (args.comparisons / 'domain_comparisons.tsv').open() as handle:
        groups = defaultdict(list)
        for row in csv.DictReader(handle, delimiter='\t'):
            groups[(row['marker'], row['taxon_a'], row['taxon_b'])].append(row)
    links = {(r['marker'], r['taxon_id']): r for r in csv.DictReader((sources['snapshot'] / 'marker_structure_links.tsv').open(), delimiter='\t')}
    needed = {key for marker, a, b in groups if len(groups[(marker, a, b)]) >= 2 for key in [(marker, a), (marker, b)]}
    mapped = defaultdict(dict)
    with gzip.open(sources['snapshot'] / 'matrix_to_structure_residues.tsv.gz', 'rt') as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            key = (r['marker'], r['taxon_id'])
            if key in needed:
                mapped[key][int(r['matrix_column_1based'])] = (int(r['protein_residue_1based']), float(r['ca_plddt']))
    matrix_dir = ROOT / 'results/phylogeny/profile-matrix-50-v1'
    checked_receipt(matrix_dir)
    mapping_receipt = json.loads((sources['snapshot'] / 'receipt.json').read_text())
    if sha(matrix_dir / 'receipt.json') != mapping_receipt['matrix_receipt_sha256']:
        raise ValueError('Changed matrix receipt')
    matrix = {r.id: str(r.seq) for r in SeqIO.parse(matrix_dir / 'matrix.faa', 'fasta')}
    coords, paes = {}, {}
    for key in needed:
        row = links[key]
        path = ROOT / row['model_path']
        if row['model_path'] not in coords:
            if sha(path) != row['model_sha256']:
                raise ValueError('Changed coordinates')
            cif = MMCIF2Dict(str(path))
            coords[row['model_path']] = {int(p): np.array([float(x), float(y), float(z)]) for atom, p, x, y, z in zip(
                cif['_atom_site.label_atom_id'], cif['_atom_site.label_seq_id'], cif['_atom_site.Cartn_x'],
                cif['_atom_site.Cartn_y'], cif['_atom_site.Cartn_z']) if atom == 'CA'}
    needed_models = {(links[k]['model_id'], links[k]['model_version']) for k in needed}
    for row in json.loads((sources['pae'] / 'pae_manifest.json').read_text()):
        key = (row['model_id'], str(row['version']))
        if row['status'] == 'verified' and key in needed_models:
            path = ROOT / row['path']
            if sha(path) != row['gzip_sha256']:
                raise ValueError('Changed PAE')
            paes[key] = (validate_pae(gzip.decompress(path.read_bytes()), row['length']), row['sequence_sha256'])
    results = []
    for (marker, a, b), domains in sorted(groups.items()):
        if len(domains) < 2:
            continue
        la, lb = links[(marker, a)], links[(marker, b)]
        ma, mb = mapped[(marker, a)], mapped[(marker, b)]
        pae_a, hash_a = paes[(la['model_id'], la['model_version'])]
        pae_b, hash_b = paes[(lb['model_id'], lb['model_version'])]
        if hash_a != la['sequence_sha256'] or hash_b != lb['sequence_sha256']:
            raise ValueError('PAE sequence identity differs')
        sites = {}
        for d in domains:
            cols = [c for c in sorted(ma.keys() & mb.keys()) if ma[c][1] >= 70 and mb[c][1] >= 70
                and matrix[a][c - 1] in AA and matrix[b][c - 1] in AA
                and int(d['domain_a_start']) <= ma[c][0] <= int(d['domain_a_end'])
                and int(d['domain_b_start']) <= mb[c][0] <= int(d['domain_b_end'])]
            if len(cols) != int(d['domain_compared_residues']):
                raise ValueError('Domain matched positions not reproduced')
            sites[d['pfam_accession']] = (np.array([ma[c][0] for c in cols]), np.array([mb[c][0] for c in cols]))
        for d1, d2 in itertools.combinations(sorted(domains, key=lambda d: d['pfam_accession']), 2):
            p1, q1 = sites[d1['pfam_accession']]
            p2, q2 = sites[d2['pfam_accession']]
            if set(p1) & set(p2) or set(q1) & set(q2):
                raise ValueError('Supposedly unambiguous domains overlap')
            xa, xb = [np.array([coords[la['model_path']][p] for p in positions]) for positions in [p1, p2]]
            ya, yb = [np.array([coords[lb['model_path']][q] for q in positions]) for positions in [q1, q2]]
            delta = np.abs(np.linalg.norm(xa[:, None] - xb[None, :], axis=2) - np.linalg.norm(ya[:, None] - yb[None, :], axis=2))
            nonadjacent = (np.abs(p1[:, None] - p2[None, :]) >= 3) & (np.abs(q1[:, None] - q2[None, :]) >= 3)
            for cutoff in [5, 10, 15]:
                confident = confidence_mask(pae_a, pae_b, np.concatenate([p1, p2]), np.concatenate([q1, q2]), cutoff)[:len(p1), len(p1):]
                good = nonadjacent & confident
                results.append({'marker': marker, 'taxon_a': a, 'taxon_b': b,
                    'domain_1': d1['pfam_accession'], 'domain_2': d2['pfam_accession'],
                    'domain_1_name': d1['pfam_name'], 'domain_2_name': d2['pfam_name'],
                    'domain_1_residues': len(p1), 'domain_2_residues': len(p2),
                    'domain_1_own_fit_rmsd': float(d1['domain_own_fit_rmsd_angstrom']),
                    'domain_2_own_fit_rmsd': float(d2['domain_own_fit_rmsd_angstrom']),
                    'domain_1_pae10_local_fraction': d1['pae10_local_pair_fraction'],
                    'domain_2_pae10_local_fraction': d2['pae10_local_pair_fraction'],
                    'pae_cutoff_angstrom': cutoff, 'interdomain_residue_pairs': int(nonadjacent.sum()),
                    'confident_interdomain_pairs': int(good.sum()),
                    'confident_interdomain_fraction': float(good.sum() / nonadjacent.sum()) if nonadjacent.any() else '',
                    'interdomain_mean_absolute_distance_change_angstrom': float(delta[nonadjacent].mean()) if nonadjacent.any() else '',
                    'confident_interdomain_mean_absolute_distance_change_angstrom': float(delta[good].mean()) if good.any() else ''})
    args.output.mkdir(parents=True)
    path = args.output / 'interdomain_confidence.tsv'
    with path.open('w') as handle:
        if results:
            writer = csv.DictWriter(handle, list(results[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(results)
    result = {'status': 'complete_interdomain_confidence_snapshot', 'rows': len(results),
        'taxon_pairs': len({(r['marker'], r['taxon_a'], r['taxon_b']) for r in results}),
        'domain_pair_taxon_combinations': len(results) // 3,
        'source_receipt_sha256': sha(args.comparisons / 'receipt.json'), 'script_sha256': sha(Path(__file__)),
        'interpretation': 'All nonadjacent matched cross-domain residue pairs, both PAE directions in both models. Low PAE-filter retention is prediction uncertainty, not proof of biological flexibility or annotation error. Differences under changing masks are descriptive and not independent evolutionary events.',
        'artifacts': {path.name: sha(path)}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
