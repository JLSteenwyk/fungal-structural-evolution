#!/usr/bin/env python3
"""Compare conserved domain geometry while retaining ambiguity and global-fit effects."""
import argparse
import csv
import gzip
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from compare_marker_structures import AA, ROOT, geometry, sha
from assess_pae_sensitivity import checked_receipt, confidence_mask
from retrieve_marker_pae import validate_pae


def fit_residuals(x, y):
    """Proper rigid fit; return per-residue squared displacement in Angstrom squared."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.shape != y.shape or x.ndim != 2 or x.shape[1] != 3 or len(x) < 3:
        raise ValueError('At least three paired CA coordinates required')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Nonfinite coordinates')
    xc, yc = x - x.mean(0), y - y.mean(0)
    u, _, vt = np.linalg.svd(xc.T @ yc)
    rotation = u @ np.diag([1., 1., 1. if np.linalg.det(u @ vt) >= 0 else -1.]) @ vt
    return np.sum((xc @ rotation - yc) ** 2, axis=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--comparisons', type=Path, required=True)
    parser.add_argument('--annotations', type=Path, required=True)
    parser.add_argument('--pae', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipts = {name: checked_receipt(getattr(args, name)) for name in ['snapshot', 'comparisons', 'annotations', 'pae']}
    for name in ['comparisons', 'pae']:
        if receipts[name]['mapping_receipt_sha256'] != sha(args.snapshot / 'receipt.json'):
            raise ValueError('Inconsistent structural snapshots')
    marker_inputs = ROOT / 'data/domains/marker-inputs-v1'
    checked_receipt(marker_inputs)
    search = ROOT / 'results/domains/marker-search-v1'
    if receipts['annotations']['search_receipt_sha256'] != sha(search / 'receipt.json'):
        raise ValueError('Changed annotation source search')
    if json.loads((search / 'config.json').read_text())['input_receipt_sha256'] != sha(marker_inputs / 'receipt.json'):
        raise ValueError('Annotations and marker sequence inputs disagree')
    if args.output.exists():
        raise FileExistsError('Use a new immutable comparison snapshot')
    with (args.snapshot / 'marker_structure_links.tsv').open() as handle:
        links = {(r['marker'], r['taxon_id']): r for r in csv.DictReader(handle, delimiter='\t')}
    with (marker_inputs / 'protein_links.tsv').open() as handle:
        protein_links = {(r['marker'], r['taxon_id']): r for r in csv.DictReader(handle, delimiter='\t')}
    for key, link in links.items():
        if any(link[k] != protein_links[key][k] for k in ['protein_id', 'sequence_sha256']):
            raise ValueError('Domain/structure full-sequence or protein identity mismatch')
    by_sequence = defaultdict(list)
    with (args.annotations / 'raw_annotated_hits.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            for field in ['alignment_start', 'alignment_end', 'protein_length']:
                row[field] = int(row[field])
            by_sequence[row['sequence_id']].append(row)
    overlapping = set()
    with (args.annotations / 'overlapping_hits.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            overlapping.update([row['hit_a'], row['hit_b']])
    eligible, domain_audit = {}, []
    for key, link in links.items():
        rows = by_sequence[protein_links[key]['sequence_id']]
        copies = Counter(r['pfam_accession'] for r in rows)
        selected = {}
        for row in rows:
            reasons = []
            if row['pfam_type'] != 'Domain':
                reasons.append('Pfam_type_not_Domain')
            if copies[row['pfam_accession']] != 1:
                reasons.append('multiple_instances_of_family')
            if row['hit_id'] in overlapping:
                reasons.append('overlaps_another_GA_hit')
            if float(row['hmm_coverage']) < .5:
                reasons.append('less_than_half_HMM_coverage')
            domain_audit.append({'marker': key[0], 'taxon_id': key[1], 'protein_id': link['protein_id'],
                'hit_id': row['hit_id'], 'pfam_accession': row['pfam_accession'],
                'pfam_type': row['pfam_type'], 'alignment_start': row['alignment_start'],
                'alignment_end': row['alignment_end'], 'status': ';'.join(reasons) if reasons else 'eligible'})
            if not reasons:
                selected[row['pfam_accession']] = row
        eligible[key] = selected
    mapped = defaultdict(dict)
    with gzip.open(args.snapshot / 'matrix_to_structure_residues.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            mapped[(row['marker'], row['taxon_id'])][int(row['matrix_column_1based'])] = (
                int(row['protein_residue_1based']), float(row['ca_plddt']))
    matrix_dir = ROOT / 'results/phylogeny/profile-matrix-50-v1'
    checked_receipt(matrix_dir)
    if sha(matrix_dir / 'receipt.json') != receipts['snapshot']['matrix_receipt_sha256']:
        raise ValueError('Changed matrix source')
    matrix = {r.id: str(r.seq) for r in SeqIO.parse(matrix_dir / 'matrix.faa', 'fasta')}
    coordinates, paes = {}, {}
    for link in links.values():
        if link['model_path'] not in coordinates:
            path = ROOT / link['model_path']
            if sha(path) != link['model_sha256']:
                raise ValueError('Changed coordinates')
            cif = MMCIF2Dict(str(path))
            coordinates[link['model_path']] = {int(p): np.array([float(x), float(y), float(z)])
                for atom, p, x, y, z in zip(cif['_atom_site.label_atom_id'], cif['_atom_site.label_seq_id'],
                cif['_atom_site.Cartn_x'], cif['_atom_site.Cartn_y'], cif['_atom_site.Cartn_z']) if atom == 'CA'}
    for row in json.loads((args.pae / 'pae_manifest.json').read_text()):
        if row['status'] == 'verified':
            path = ROOT / row['path']
            if sha(path) != row['gzip_sha256']:
                raise ValueError('Changed PAE')
            paes[(row['model_id'], str(row['version']))] = (validate_pae(gzip.decompress(path.read_bytes()), row['length']), row)
    with (args.comparisons / 'pairwise_metrics.tsv').open() as handle:
        pairs = [r for r in csv.DictReader(handle, delimiter='\t') if r['plddt_cutoff'] == '70']
    results, exclusions = [], []
    for baseline in pairs:
        ka, kb = (baseline['marker'], baseline['taxon_a']), (baseline['marker'], baseline['taxon_b'])
        a, b = links[ka], links[kb]
        ma, mb = mapped[ka], mapped[kb]
        common = sorted(set(ma) & set(mb))
        cols = [c for c in common if ma[c][1] >= 70 and mb[c][1] >= 70
                and matrix[ka[1]][c - 1] in AA and matrix[kb[1]][c - 1] in AA]
        px, py = np.array([ma[c][0] for c in cols]), np.array([mb[c][0] for c in cols])
        x = np.array([coordinates[a['model_path']][p] for p in px])
        y = np.array([coordinates[b['model_path']][p] for p in py])
        global_squared = fit_residuals(x, y)
        if len(cols) != int(baseline['compared_residues']) or not np.isclose(np.sqrt(global_squared.mean()), float(baseline['ca_superposition_rmsd_angstrom']), atol=1e-8):
            raise ValueError('Original whole-marker comparison not reproduced')
        base = {k: baseline[k] for k in ['marker', 'taxon_a', 'taxon_b', 'same_model_coordinates']}
        families = sorted(eligible[ka].keys() & eligible[kb].keys())
        if not families:
            exclusions.append(base | {'pfam_accession': '', 'reason': 'no_shared_eligible_single_instance_domain'})
        for family in families:
            da, db = eligible[ka][family], eligible[kb][family]
            inside = lambda c: (da['alignment_start'] <= ma[c][0] <= da['alignment_end']
                                and db['alignment_start'] <= mb[c][0] <= db['alignment_end'])
            domain_common = [c for c in common if inside(c)]
            indices = np.array([i for i, c in enumerate(cols) if inside(c)], dtype=int)
            if len(indices) < 30 or len(indices) < .5 * len(domain_common):
                exclusions.append(base | {'pfam_accession': family, 'reason': 'fewer_than_30_qualified_or_half_domain_shared_positions'})
                continue
            dx, dy, pa, pb = x[indices], y[indices], px[indices], py[indices]
            g = geometry(dx, dy, pa, pb)
            fitted_global = float(np.sqrt(global_squared[indices].mean()))
            fitted_domain = g['ca_superposition_rmsd_angstrom']
            if fitted_domain > fitted_global + 1e-7:
                raise ValueError('Domain optimal fit is worse than the same sites under a global fit')
            p1, p2 = paes[(a['model_id'], a['model_version'])], paes[(b['model_id'], b['model_version'])]
            if p1[1]['sequence_sha256'] != a['sequence_sha256'] or p2[1]['sequence_sha256'] != b['sequence_sha256']:
                raise ValueError('PAE sequence identity differs')
            confident = confidence_mask(p1[0], p2[0], pa, pb, 10)
            dist_x = np.linalg.norm(dx[:, None] - dx[None, :], axis=2)
            dist_y = np.linalg.norm(dy[:, None] - dy[None, :], axis=2)
            local = (np.triu(np.ones(dist_x.shape, bool), 1) & ((dist_x <= 15) | (dist_y <= 15))
                     & (np.abs(pa[:, None] - pa[None, :]) >= 3) & (np.abs(pb[:, None] - pb[None, :]) >= 3))
            good = local & confident
            if int(local.sum()) != g['local_distance_pairs']:
                raise ValueError('Local residue pairs differ')
            results.append(base | {'pfam_accession': family, 'pfam_name': da['pfam_name'],
                'hit_a': da['hit_id'], 'hit_b': db['hit_id'],
                'domain_a_start': da['alignment_start'], 'domain_a_end': da['alignment_end'],
                'domain_b_start': db['alignment_start'], 'domain_b_end': db['alignment_end'],
                'domain_shared_positions': len(domain_common), 'domain_compared_residues': len(indices),
                'whole_marker_compared_residues': len(cols),
                'whole_marker_rmsd_angstrom': float(baseline['ca_superposition_rmsd_angstrom']),
                'domain_sites_under_whole_marker_fit_rmsd_angstrom': fitted_global,
                'domain_own_fit_rmsd_angstrom': fitted_domain,
                'domain_fit_improvement_angstrom': fitted_global - fitted_domain,
                'domain_uncorrected_sequence_difference': float(np.mean([matrix[ka[1]][cols[i] - 1] != matrix[kb[1]][cols[i] - 1] for i in indices])),
                'local_distance_pairs': int(local.sum()),
                'local_mean_absolute_distance_change_angstrom': g['local_distance_mean_absolute_change_angstrom'],
                'pae10_local_pairs': int(good.sum()),
                'pae10_local_pair_fraction': float(good.sum() / local.sum()) if local.any() else '',
                'pae10_local_mean_absolute_distance_change_angstrom': float(np.abs(dist_x - dist_y)[good].mean()) if good.any() else ''})
    args.output.mkdir(parents=True)
    for name, rows in [('domain_comparisons.tsv', results), ('excluded_pairs.tsv', exclusions), ('domain_eligibility.tsv', domain_audit)]:
        with (args.output / name).open('w') as handle:
            if rows:
                writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
                writer.writeheader(); writer.writerows(rows)
    receipt = {'status': 'complete_domain_comparison_snapshot',
        'sources': {k: {'path': str(getattr(args, k)), 'receipt_sha256': sha(getattr(args, k) / 'receipt.json')} for k in receipts},
        'script_sha256': sha(Path(__file__)), 'geometry_script_sha256': sha(ROOT / 'scripts/compare_marker_structures.py'),
        'full_marker_pairs_screened': len(pairs), 'domain_comparisons': len(results),
        'distinct_marker_taxon_pairs': len({(r['marker'], r['taxon_a'], r['taxon_b']) for r in results}),
        'distinct_pfam_domains': len({r['pfam_accession'] for r in results}),
        'excluded_pair_or_domain_rows': len(exclusions), 'mapped_annotation_hits': len(domain_audit),
        'eligible_annotation_hits': sum(r['status'] == 'eligible' for r in domain_audit),
        'filters': 'Pfam Domain type; single instance; no overlap with any GA hit; >=50% HMM coverage; matched profile sites in both domain spans; both pLDDT>=70; >=30 qualified residues and >=50% shared domain sites; same direct-comparison baseline.',
        'interpretation': 'Descriptive conserved-domain geometry. Separate fitting necessarily improves fit and is not proof of biological domain motion. Pairwise rows are dependent; no branch rates, gains/losses, selection or significance are inferred. Restricted to existing immutable AFDB comparison snapshot.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
