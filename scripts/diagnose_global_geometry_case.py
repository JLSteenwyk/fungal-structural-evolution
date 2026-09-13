#!/usr/bin/env python3
"""Inspect the largest RMSD benchmark pair using shared Pfam regions and PAE."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from assess_pae_sensitivity import checked_receipt, confidence_mask
from compare_marker_structures import ROOT, sha
from compare_marker_domains import fit_residuals
from prepare_paired_phylogenetic_inputs import write_table
from retrieve_marker_pae import validate_pae


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['benchmark', 'inputs', 'snapshot', 'annotations', 'pae', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    receipts = {name: checked_receipt(getattr(args, name)) for name in ['benchmark', 'inputs', 'snapshot', 'annotations', 'pae']}
    for name in ['inputs', 'snapshot', 'pae']:
        if receipts['benchmark']['source_receipts'][name]['sha256'] != sha(getattr(args, name) / 'receipt.json'):
            raise ValueError('Case sources differ from benchmark')
    if args.output.exists():
        raise FileExistsError('Use a new immutable case output')
    rows = list(csv.DictReader((args.benchmark / 'path_geometry.tsv').open(), delimiter='\t'))
    focal = max(rows, key=lambda r: float(r['ca_superposition_rmsd_angstrom']))
    marker, a, b = focal['marker'], focal['taxon_a'], focal['taxon_b']
    taxa = [a, b]
    links = {r['taxon_id']: r for r in csv.DictReader((args.snapshot / 'marker_structure_links.tsv').open(), delimiter='\t')
             if r['marker'] == marker and r['taxon_id'] in taxa}
    hits = defaultdict(list)
    for row in csv.DictReader((args.annotations / 'raw_annotated_hits.tsv').open(), delimiter='\t'):
        for taxon in taxa:
            if row['sequence_id'] == 'S' + links[taxon]['sequence_sha256']:
                hits[taxon].append(row)
    overlap = {r[k] for r in csv.DictReader((args.annotations / 'overlapping_hits.tsv').open(), delimiter='\t') for k in ['hit_a', 'hit_b']}
    eligible = {}
    for taxon in taxa:
        copies = Counter(r['pfam_accession'] for r in hits[taxon])
        eligible[taxon] = {r['pfam_accession']: r for r in hits[taxon] if copies[r['pfam_accession']] == 1
                          and r['hit_id'] not in overlap and float(r['hmm_coverage']) >= .5}
    sequences = {r.id: str(r.seq) for r in SeqIO.parse(args.inputs / marker / 'aa.faa', 'fasta')}
    source_columns = [int(r['matrix_column_1based']) for r in csv.DictReader((args.inputs / marker / 'columns.tsv').open(), delimiter='\t')]
    common = [i for i in range(len(source_columns)) if sequences[a][i] != '?' and sequences[b][i] != '?']
    if len(common) != int(focal['geometry_matched_residues']):
        raise ValueError('Case sites differ from benchmark')
    mapped = defaultdict(dict)
    with gzip.open(args.snapshot / 'matrix_to_structure_residues.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['marker'] == marker and row['taxon_id'] in taxa:
                mapped[row['taxon_id']][int(row['matrix_column_1based'])] = int(row['protein_residue_1based'])
    pae_manifest = {(r['model_id'], str(r['version'])): r for r in json.loads((args.pae / 'pae_manifest.json').read_text())}
    positions, coordinates, errors = {}, {}, {}
    for taxon in taxa:
        link = links[taxon]
        path = ROOT / link['model_path']
        if sha(path) != link['model_sha256']:
            raise ValueError('Changed case model')
        cif = MMCIF2Dict(str(path))
        ca = {int(p): [float(x), float(y), float(z)] for atom, p, x, y, z in zip(cif['_atom_site.label_atom_id'],
              cif['_atom_site.label_seq_id'], cif['_atom_site.Cartn_x'], cif['_atom_site.Cartn_y'], cif['_atom_site.Cartn_z']) if atom == 'CA'}
        positions[taxon] = np.array([mapped[taxon][source_columns[i]] for i in common])
        coordinates[taxon] = np.array([ca[p] for p in positions[taxon]])
        entry = pae_manifest[link['model_id'], link['model_version']]
        path = ROOT / entry['path']
        if sha(path) != entry['gzip_sha256'] or entry['sequence_sha256'] != link['sequence_sha256']:
            raise ValueError('Changed case PAE source')
        raw = gzip.decompress(path.read_bytes())
        if hashlib.sha256(raw).hexdigest() != entry['json_sha256']:
            raise ValueError('Changed PAE data')
        errors[taxon] = validate_pae(raw, entry['length'])
    x, y, pa, pb = coordinates[a], coordinates[b], positions[a], positions[b]
    residuals = fit_residuals(x, y)
    if not np.isclose(np.sqrt(residuals.mean()), float(focal['ca_superposition_rmsd_angstrom']), atol=1e-8):
        raise ValueError('Global case fit differs from benchmark')
    dx = np.linalg.norm(x[:, None] - x[None, :], axis=2)
    dy = np.linalg.norm(y[:, None] - y[None, :], axis=2)
    pae_mask = confidence_mask(errors[a], errors[b], pa, pb, 10)
    region_rows = []
    for family in sorted(set(eligible[a]) & set(eligible[b])):
        ha, hb = eligible[a][family], eligible[b][family]
        inside = ((pa >= int(ha['alignment_start'])) & (pa <= int(ha['alignment_end']))
                  & (pb >= int(hb['alignment_start'])) & (pb <= int(hb['alignment_end'])))
        if inside.sum() < 30:
            continue
        outside = ((pa < int(ha['alignment_start'])) | (pa > int(ha['alignment_end']))) & ((pb < int(hb['alignment_start'])) | (pb > int(hb['alignment_end'])))
        cross = (inside[:, None] & outside[None, :] & (np.abs(pa[:, None] - pa[None, :]) >= 3)
                 & (np.abs(pb[:, None] - pb[None, :]) >= 3))
        confident = cross & pae_mask
        within = np.triu(inside[:, None] & inside[None, :], 1)
        row = {'marker': marker, 'taxon_a': a, 'taxon_b': b, 'pfam_accession': family,
            'pfam_name': ha['pfam_name'], 'pfam_type': ha['pfam_type'],
            'matched_region_residues': int(inside.sum()), 'global_marker_rmsd_angstrom': float(np.sqrt(residuals.mean())),
            'region_under_global_fit_rmsd_angstrom': float(np.sqrt(residuals[inside].mean())),
            'region_own_fit_rmsd_angstrom': float(np.sqrt(fit_residuals(x[inside], y[inside]).mean())),
            'within_region_pair_pae10_fraction': float((within & pae_mask).sum() / within.sum()),
            'cross_region_pairs': int(cross.sum()), 'cross_region_pae10_pairs': int(confident.sum()),
            'cross_region_pae10_fraction': float(confident.sum() / cross.sum()) if cross.any() else '',
            'cross_region_unfiltered_mean_distance_change_angstrom': float(np.mean(np.abs(dx[cross] - dy[cross]))) if cross.any() else '',
            'cross_region_pae10_mean_distance_change_angstrom': float(np.mean(np.abs(dx[confident] - dy[confident]))) if confident.any() else ''}
        region_rows.append(row)
    args.output.mkdir(parents=True)
    write_table(args.output / 'region_geometry.tsv', region_rows)
    write_table(args.output / 'pfam_hits.tsv', [dict(taxon_id=t, **r) for t in taxa for r in hits[t]])
    result = {'status': 'complete_maximum_rmsd_case_diagnostic', 'selection': 'Largest whole-marker RMSD among all accepted tree-path benchmark pairs',
        'marker': marker, 'taxon_a': a, 'taxon_b': b, 'benchmark_row': focal,
        'regions_compared': len(region_rows),
        'source_receipts': {name: {'path': str(getattr(args, name)), 'sha256': sha(getattr(args, name) / 'receipt.json')} for name in receipts},
        'script_sha256': sha(Path(__file__)),
        'interpretation': 'Single-instance, nonoverlapping, >=half-HMM shared Pfam regions are compared on identical retained sites. Repeat/Family annotations are not reclassified as independent domains. An improved separate rigid fit is mathematically expected and alone does not prove biological motion. Cross-region PAE qualifies orientation interpretability. Missing GA hits do not prove domain loss, and hybrid/taxonomic review remains necessary.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(region_rows, indent=2))


if __name__ == '__main__':
    main()
