#!/usr/bin/env python3
"""Benchmark fitted tree path distances against matched-residue CA geometry."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import numpy as np
from Bio import Phylo, SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from assess_pae_sensitivity import checked_receipt, confidence_mask
from compare_marker_structures import ROOT, sha, geometry
from prepare_paired_phylogenetic_inputs import write_table
from retrieve_marker_pae import validate_pae
from run_paired_marker_fits import tree_edges


def path_mask(splits, a, b):
    if a == b:
        raise ValueError('Distinct tips required')
    return np.array([(a in split) != (b in split) for split in splits])


def path_interval(branch_draws, mask):
    draws = np.asarray(branch_draws, dtype=float)
    if draws.ndim != 2 or draws.shape[1] != len(mask) or not np.isfinite(draws).all() or (draws < 0).any():
        raise ValueError('Invalid joint branch draws')
    # Sum within each jointly fitted replicate, not marginal interval endpoints.
    return np.quantile(draws[:, mask].sum(axis=1), [.025, .5, .975], method='linear')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['inputs', 'snapshot', 'pae', 'fits', 'fit_audit', 'resampling', 'resampling_audit', 'output']:
        parser.add_argument('--' + name.replace('_', '-'), type=Path, required=True)
    args = parser.parse_args()
    checked = {name: checked_receipt(getattr(args, name)) for name in ['inputs', 'snapshot', 'pae', 'fit_audit', 'resampling_audit']}
    if (checked['inputs']['source_receipts']['snapshot']['sha256'] != sha(args.snapshot / 'receipt.json')
            or checked['pae']['mapping_receipt_sha256'] != sha(args.snapshot / 'receipt.json')
            or checked['fit_audit']['fit_receipt_sha256'] != sha(args.fits / 'receipt.json')
            or checked['resampling_audit']['resampling_receipt_sha256'] != sha(args.resampling / 'receipt.json')):
        raise ValueError('Mismatched source snapshots')
    fit_config = json.loads((args.fits / 'config.json').read_text())
    resampling_config = json.loads((args.resampling / 'config.json').read_text())
    fit_receipt = json.loads((args.fits / 'receipt.json').read_text())
    resampling_receipt = json.loads((args.resampling / 'receipt.json').read_text())
    if (fit_config['input_receipt_sha256'] != sha(args.inputs / 'receipt.json')
            or fit_receipt['config_sha256'] != sha(args.fits / 'config.json')
            or resampling_receipt['config_sha256'] != sha(args.resampling / 'config.json')
            or resampling_config['source_receipts']['fits'] != sha(args.fits / 'receipt.json')
            or resampling_config['source_receipts']['inputs'] != sha(args.inputs / 'receipt.json')):
        raise ValueError('Fit input lineage differs')
    if args.output.exists():
        raise FileExistsError('Use a new immutable benchmark output')
    links = {(r['marker'], r['taxon_id']): r for r in csv.DictReader((args.snapshot / 'marker_structure_links.tsv').open(), delimiter='\t')}
    mapped = defaultdict(dict)
    with gzip.open(args.snapshot / 'matrix_to_structure_residues.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            mapped[row['marker'], row['taxon_id']][int(row['matrix_column_1based'])] = int(row['protein_residue_1based'])
    pae_rows = {(r['model_id'], str(r['version'])): r for r in json.loads((args.pae / 'pae_manifest.json').read_text())}
    coordinates, errors = {}, {}
    def load_model(link):
        name = link['model_path']
        if name not in coordinates:
            path = ROOT / name
            if sha(path) != link['model_sha256']:
                raise ValueError('Changed coordinates')
            cif = MMCIF2Dict(str(path))
            ca = {}
            for atom, pos, x, y, z in zip(cif['_atom_site.label_atom_id'], cif['_atom_site.label_seq_id'],
                    cif['_atom_site.Cartn_x'], cif['_atom_site.Cartn_y'], cif['_atom_site.Cartn_z']):
                if atom == 'CA':
                    if int(pos) in ca:
                        raise ValueError('Ambiguous CA atom')
                    ca[int(pos)] = np.array([float(x), float(y), float(z)])
            coordinates[name] = ca
            row = pae_rows[link['model_id'], link['model_version']]
            if row['sequence_sha256'] != link['sequence_sha256']:
                raise ValueError('PAE sequence identity differs')
            compressed = (ROOT / row['path']).read_bytes()
            if hashlib.sha256(compressed).hexdigest() != row['gzip_sha256']:
                raise ValueError('Changed PAE file')
            raw = gzip.decompress(compressed)
            if hashlib.sha256(raw).hexdigest() != row['json_sha256']:
                raise ValueError('Changed PAE content')
            errors[name] = validate_pae(raw, row['length'])
        return coordinates[name], errors[name]
    batches = {(r['marker'], r['block_length']): r for r in resampling_receipt['results']}
    results, exclusions = [], []
    for fitted_marker in fit_receipt['results']:
        marker = fitted_marker['marker']
        folder = args.inputs / marker
        aa = {r.id: str(r.seq) for r in SeqIO.parse(folder / 'aa.faa', 'fasta')}
        states = {r.id: str(r.seq) for r in SeqIO.parse(folder / '3di.faa', 'fasta')}
        if set(aa) != set(states):
            raise ValueError('Paired taxa differ')
        for taxon in aa:
            if [c == '?' for c in aa[taxon]] != [c == '?' for c in states[taxon]]:
                raise ValueError('Paired missingness differs')
        source_columns = [int(r['matrix_column_1based']) for r in csv.DictReader((folder / 'columns.tsv').open(), delimiter='\t')]
        trees, edges = {}, {}
        for label in ['aa', '3di_af', '3di_af_empirical', '3di_llm']:
            path = args.fits / marker / (label + '.treefile')
            r = json.loads((args.fits / marker / (label + '.receipt.json')).read_text())
            if sha(path) != r['artifacts'][label + '.treefile']:
                raise ValueError('Changed fitted tree')
            trees[label] = Phylo.read(path, 'newick')
            edges[label] = tree_edges(path, set(aa))
        splits = sorted(edges['aa'])
        if any(set(e) != set(splits) for e in edges.values()):
            raise ValueError('Fitted topologies differ')
        batch = batches[marker, 30]
        batch_folder = args.resampling / marker / 'block-30'
        if sha(batch_folder / 'receipt.json') != batch['receipt_sha256']:
            raise ValueError('Changed joint-resampling batch')
        br = json.loads((batch_folder / 'receipt.json').read_text())
        draws = {'aa': [], '3di': []}
        for name, checksum in br['receipt_hashes'].items():
            path = batch_folder / name
            if sha(path) != checksum:
                raise ValueError('Changed audited joint branch draws')
            draw = json.loads(path.read_text())
            if draw['status'] == 'unestimable_draw':
                continue
            if draw['status'] != 'completed_draw' or draw['splits'] != [list(s) for s in splits]:
                raise ValueError('Unexpected resampling state or splits')
            for label in draws:
                draws[label].append(draw['fits'][label]['branch_lengths'])
        if len(draws['aa']) != br['completed_draws'] or len(draws['aa']) < .9 * resampling_config['replicates']:
            raise ValueError('Insufficient audited joint draws')
        for a, b in combinations(sorted(aa), 2):
            common = [i for i in range(len(source_columns)) if aa[a][i] != '?' and aa[b][i] != '?']
            observed_a = sum(c != '?' for c in aa[a]); observed_b = sum(c != '?' for c in aa[b])
            row = {'marker': marker, 'taxon_a': a, 'taxon_b': b, 'tree_taxa': len(aa),
                'tree_alignment_columns': len(source_columns), 'observed_a': observed_a, 'observed_b': observed_b,
                'geometry_matched_residues': len(common), 'fraction_of_tree_columns': len(common) / len(source_columns),
                'same_model_coordinates': links[marker, a]['model_sha256'] == links[marker, b]['model_sha256']}
            mask = path_mask(splits, a, b)
            for label in edges:
                distance = sum(edges[label][split] for split, keep in zip(splits, mask) if keep)
                if not np.isclose(distance, trees[label].distance(a, b), rtol=0, atol=1e-10):
                    raise ValueError('Split-based and tree-traversal path distances differ')
                row[label + '_path_distance'] = distance
            for label in draws:
                q = path_interval(draws[label], mask)
                for suffix, value in zip(['p2_5', 'median', 'p97_5'], q):
                    row[label + '_path_block30_' + suffix] = float(value)
            row['estimable_block30_draws'] = len(draws['aa'])
            if len(common) < 50 or len(common) < .5 * max(observed_a, observed_b):
                exclusions.append(row | {'reason': 'fewer_than_50_shared_sites_or_below_half_of_either_taxon_observations'})
                continue
            la, lb = links[marker, a], links[marker, b]
            ca, pae_a = load_model(la); cb, pae_b = load_model(lb)
            pa = np.array([mapped[marker, a][source_columns[i]] for i in common])
            pb = np.array([mapped[marker, b][source_columns[i]] for i in common])
            x = np.array([ca[p] for p in pa]); y = np.array([cb[p] for p in pb])
            row.update(geometry(x, y, pa, pb))
            row['uncorrected_aa_difference'] = sum(aa[a][i] != aa[b][i] for i in common) / len(common)
            row['uncorrected_3di_difference'] = sum(states[a][i] != states[b][i] for i in common) / len(common)
            dx = np.linalg.norm(x[:, None] - x[None, :], axis=2)
            dy = np.linalg.norm(y[:, None] - y[None, :], axis=2)
            local = (np.triu(np.ones(dx.shape, dtype=bool), 1) & ((dx <= 15) | (dy <= 15))
                     & (np.abs(pa[:, None] - pa[None, :]) >= 3) & (np.abs(pb[:, None] - pb[None, :]) >= 3))
            retained = local & confidence_mask(pae_a, pae_b, pa, pb, 10)
            row['pae10_local_pairs'] = int(retained.sum())
            row['pae10_local_retained_fraction'] = float(retained.sum() / local.sum()) if local.any() else ''
            row['pae10_local_mean_absolute_change_angstrom'] = float(np.mean(np.abs(dx[retained] - dy[retained]))) if retained.any() else ''
            results.append(row)
        print(marker, 'benchmarked', flush=True)
    args.output.mkdir(parents=True)
    write_table(args.output / 'path_geometry.tsv', results)
    if exclusions:
        write_table(args.output / 'coverage_exclusions.tsv', exclusions)
    receipt = {'status': 'complete_tree_path_geometry_benchmark', 'markers': len(fit_receipt['results']),
        'accepted_taxon_pairs': len(results), 'excluded_taxon_pairs': len(exclusions),
        'source_receipts': {name: {'path': str(getattr(args, name)), 'sha256': sha(getattr(args, name) / 'receipt.json')}
                            for name in ['inputs', 'snapshot', 'pae', 'fits', 'fit_audit', 'resampling', 'resampling_audit']},
        'script_sha256': sha(Path(__file__)), 'geometry_helper_sha256': sha(Path(__file__).with_name('compare_marker_structures.py')),
        'interpretation': 'Leaf-to-leaf sums of fitted branch lengths benchmarked against matched-residue direct geometry. Path intervals sum branches within each joint block-30 draw before taking percentiles. Tree fits use all retained marker columns/taxa; pair geometry uses the jointly observed subset and records coverage. This is not a reconstruction of physical displacement on individual branches, a proof of geometric additivity, or an acceleration/coupling test. Distances and pairs share ancestry, sites and prediction sources.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'source_receipts'}, indent=2))


if __name__ == '__main__':
    main()
