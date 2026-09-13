#!/usr/bin/env python3
"""Benchmark 3Di mismatch against geometry and sequence on identical aligned residues."""
import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from compare_marker_structures import ROOT, AA, sha, geometry
from assess_pae_sensitivity import checked_receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['encodings', 'snapshot', 'comparisons', 'domains', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    receipts = {name: checked_receipt(getattr(args, name)) for name in ['encodings', 'snapshot', 'comparisons', 'domains']}
    mapping_hash = sha(args.snapshot / 'receipt.json')
    if (receipts['encodings']['mapping_receipt_sha256'] != mapping_hash
            or receipts['comparisons']['mapping_receipt_sha256'] != mapping_hash
            or receipts['domains']['sources']['snapshot']['receipt_sha256'] != mapping_hash):
        raise ValueError('Encoding and geometric benchmarks use different snapshots')
    if args.output.exists():
        raise FileExistsError('Use a new immutable benchmark output')
    encoded = {}
    for row in csv.DictReader((args.encodings / 'model_summary.tsv').open(), delimiter='\t'):
        path = ROOT / row['encoding_path']
        if sha(path) != row['encoding_sha256']:
            raise ValueError('Changed audited encoding')
        with np.load(path, allow_pickle=False) as source:
            encoded[row['model_name']] = {k: source[k].copy() for k in source.files}
        data = encoded[row['model_name']]
        data['sequence'], data['states'] = str(data['sequence']), str(data['states'])
    links = {(r['marker'], r['taxon_id']): r for r in csv.DictReader((args.snapshot / 'marker_structure_links.tsv').open(), delimiter='\t')}
    mapped = defaultdict(dict)
    with gzip.open(args.snapshot / 'matrix_to_structure_residues.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            mapped[(row['marker'], row['taxon_id'])][int(row['matrix_column_1based'])] = (int(row['protein_residue_1based']), float(row['ca_plddt']))
    matrix_dir = ROOT / 'results/phylogeny/profile-matrix-50-v1'
    checked_receipt(matrix_dir)
    if sha(matrix_dir / 'receipt.json') != receipts['snapshot']['matrix_receipt_sha256']:
        raise ValueError('Changed marker matrix')
    matrix = {r.id: str(r.seq) for r in SeqIO.parse(matrix_dir / 'matrix.faa', 'fasta')}
    coordinates = {}
    for row in links.values():
        name = Path(row['model_path']).stem
        if name in coordinates:
            continue
        path = ROOT / row['model_path']
        if sha(path) != row['model_sha256']:
            raise ValueError('Changed coordinates')
        cif = MMCIF2Dict(str(path))
        coordinates[name] = {int(p): np.array([float(x), float(y), float(z)]) for atom, p, x, y, z in zip(
            cif['_atom_site.label_atom_id'], cif['_atom_site.label_seq_id'], cif['_atom_site.Cartn_x'], cif['_atom_site.Cartn_y'], cif['_atom_site.Cartn_z']) if atom == 'CA'}
    pairs = [r for r in csv.DictReader((args.comparisons / 'pairwise_metrics.tsv').open(), delimiter='\t') if r['plddt_cutoff'] == '70']
    domains = defaultdict(list)
    for row in csv.DictReader((args.domains / 'domain_comparisons.tsv').open(), delimiter='\t'):
        domains[(row['marker'], row['taxon_a'], row['taxon_b'])].append(row)
    rows, exclusions = [], []
    regimes = ['valid_focal70', 'six_residue70', 'six_residue70_pae10']
    for pair in pairs:
        key = (pair['marker'], pair['taxon_a'], pair['taxon_b'])
        a, b = links[(key[0], key[1])], links[(key[0], key[2])]
        ma, mb = mapped[(key[0], key[1])], mapped[(key[0], key[2])]
        name_a, name_b = Path(a['model_path']).stem, Path(b['model_path']).stem
        ea, eb = encoded[name_a], encoded[name_b]
        columns = [c for c in sorted(ma.keys() & mb.keys()) if ma[c][1] >= 70 and mb[c][1] >= 70
                   and matrix[key[1]][c - 1] in AA and matrix[key[2]][c - 1] in AA]
        if len(columns) != int(pair['compared_residues']):
            raise ValueError('Whole-marker site coverage differs')
        tasks = [('whole_marker', '', columns, 50, float(pair['ca_superposition_rmsd_angstrom']))]
        for domain in domains[key]:
            dc = [c for c in columns if int(domain['domain_a_start']) <= ma[c][0] <= int(domain['domain_a_end'])
                  and int(domain['domain_b_start']) <= mb[c][0] <= int(domain['domain_b_end'])]
            if len(dc) != int(domain['domain_compared_residues']):
                raise ValueError('Domain site coverage differs')
            tasks.append(('domain', domain['pfam_accession'], dc, 30, float(domain['domain_own_fit_rmsd_angstrom'])))
        for scope, family, cols, minimum, expected_rmsd in tasks:
            pa, pb = np.array([ma[c][0] for c in cols]), np.array([mb[c][0] for c in cols])
            x = np.array([coordinates[name_a][p] for p in pa]); y = np.array([coordinates[name_b][p] for p in pb])
            baseline = geometry(x, y, pa, pb)
            if not np.isclose(baseline['ca_superposition_rmsd_angstrom'], expected_rmsd, atol=1e-8):
                raise ValueError('Baseline geometry differs')
            sa = np.array([ea['sequence'][p - 1] for p in pa]); sb = np.array([eb['sequence'][p - 1] for p in pb])
            if ''.join(sa) != ''.join(matrix[key[1]][c - 1] for c in cols) or ''.join(sb) != ''.join(matrix[key[2]][c - 1] for c in cols):
                raise ValueError('Alphabet sequence-to-alignment correspondence differs')
            if not np.allclose(ea['ca_plddt'][pa - 1], [ma[c][1] for c in cols]) or not np.allclose(eb['ca_plddt'][pb - 1], [mb[c][1] for c in cols]):
                raise ValueError('Confidence differs from structural mapping')
            ta = np.array([ea['states'][p - 1] for p in pa]); tb = np.array([eb['states'][p - 1] for p in pb])
            mask = ea['valid'][pa - 1] & eb['valid'][pb - 1]
            for regime in regimes:
                if regime == 'six_residue70':
                    mask = mask & (ea['feature_min_plddt'][pa - 1] >= 70) & (eb['feature_min_plddt'][pb - 1] >= 70)
                elif regime == 'six_residue70_pae10':
                    mask = mask & (ea['feature_max_pae'][pa - 1] <= 10) & (eb['feature_max_pae'][pb - 1] <= 10)
                count = int(mask.sum())
                base = {'marker': key[0], 'taxon_a': key[1], 'taxon_b': key[2], 'scope': scope,
                    'pfam_accession': family, 'confidence_regime': regime,
                    'original_geometric_residues': len(cols), 'compared_residues': count,
                    'retained_fraction': count / len(cols), 'same_model_coordinates': a['model_sha256'] == b['model_sha256']}
                if count < minimum or count < .5 * len(cols):
                    exclusions.append(base | {'reason': 'insufficient_qualified_3di_sites'})
                    continue
                rows.append(base | {'uncorrected_3di_difference': float(np.mean(ta[mask] != tb[mask])),
                    'uncorrected_sequence_difference': float(np.mean(sa[mask] != sb[mask])),
                    **geometry(x[mask], y[mask], pa[mask], pb[mask])})
    args.output.mkdir(parents=True)
    for name, records in [('benchmark.tsv', rows), ('exclusions.tsv', exclusions)]:
        with (args.output / name).open('w') as handle:
            if records:
                writer = csv.DictWriter(handle, list(records[0]), delimiter='\t', lineterminator='\n')
                writer.writeheader(); writer.writerows(records)
    result = {'status': 'complete_matched_site_alphabet_geometry_benchmark',
        'source_receipts': {k: {'path': str(getattr(args, k)), 'sha256': sha(getattr(args, k) / 'receipt.json')} for k in receipts},
        'script_sha256': sha(Path(__file__)), 'geometry_source_sha256': sha(ROOT / 'scripts/compare_marker_structures.py'),
        'whole_marker_pairs': len(pairs), 'domain_comparisons': sum(len(v) for v in domains.values()),
        'benchmark_rows': len(rows), 'excluded_rows': len(exclusions), 'regimes': regimes,
        'interpretation': 'Same aligned residue sets for 3Di mismatch, amino-acid mismatch and recomputed geometry. Regimes change site composition; missing or invalid states never count as similarities. 3Di distances are uncorrected state fractions, not branch lengths or Angstrom. Pairwise observations share ancestry; branch-estimate benchmarking, substitution-model assessment and source circularity controls remain pending.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
