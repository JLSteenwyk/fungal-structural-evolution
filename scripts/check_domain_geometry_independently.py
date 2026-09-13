#!/usr/bin/env python3
"""Independently recalculate one identity-hash-selected comparison per Pfam domain."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.spatial.transform import Rotation
from scipy.spatial.distance import pdist
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--domains', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use immutable output')
    receipt = checked_receipt(a.domains)
    chosen = {}
    for row in read_table(a.domains / 'domain_comparisons.tsv'):
        identity = '|'.join(row[k] for k in ['marker', 'taxon_a', 'taxon_b', 'pfam_accession'])
        digest = hashlib.sha256(('domain_geometry_audit_v1|' + identity).encode()).hexdigest()
        family = row['pfam_accession']
        if family not in chosen or digest < chosen[family][0]:
            chosen[family] = digest, row
    snapshot = Path(receipt['sources']['snapshot']['path']); sr = checked_receipt(snapshot)
    if sha(snapshot / 'receipt.json') != receipt['sources']['snapshot']['receipt_sha256']:
        raise ValueError('Snapshot changed')
    matrix_dir = Path('results/phylogeny/profile-matrix-50-v1'); checked_receipt(matrix_dir)
    if sha(matrix_dir / 'receipt.json') != sr['matrix_receipt_sha256']:
        raise ValueError('Matrix changed')
    sequences = {r.id: str(r.seq) for r in SeqIO.parse(matrix_dir / 'matrix.faa', 'fasta')}
    links = {(r['marker'], r['taxon_id']): r for r in read_table(snapshot / 'marker_structure_links.tsv')}
    needed = {(r['marker'], r[t]) for _, r in chosen.values() for t in ['taxon_a', 'taxon_b']}
    mapped = defaultdict(dict)
    with gzip.open(snapshot / 'matrix_to_structure_residues.tsv.gz', 'rt') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            key = row['marker'], row['taxon_id']
            if key in needed:
                mapped[key][int(row['matrix_column_1based'])] = int(row['protein_residue_1based']), float(row['ca_plddt'])
    coords = {}
    for key in needed:
        row = links[key]; path = Path(row['model_path'])
        if str(path) in coords:
            continue
        if sha(path) != row['model_sha256']:
            raise ValueError('Coordinate source changed')
        cif = MMCIF2Dict(str(path)); points = {}
        for atom, position, x, y, z in zip(cif['_atom_site.label_atom_id'], cif['_atom_site.label_seq_id'], cif['_atom_site.Cartn_x'], cif['_atom_site.Cartn_y'], cif['_atom_site.Cartn_z']):
            if atom == 'CA':
                if int(position) in points:
                    raise ValueError('Duplicate CA')
                points[int(position)] = [float(x), float(y), float(z)]
        coords[str(path)] = points
    records = []
    canonical = set('ACDEFGHIKLMNPQRSTVWY')
    for digest, row in chosen.values():
        ka, kb = (row['marker'], row['taxon_a']), (row['marker'], row['taxon_b'])
        ma, mb = mapped[ka], mapped[kb]
        columns = [c for c in sorted(ma.keys() & mb.keys()) if min(ma[c][1], mb[c][1]) >= 70 and sequences[ka[1]][c-1] in canonical and sequences[kb[1]][c-1] in canonical]
        indices = [i for i, c in enumerate(columns) if int(row['domain_a_start']) <= ma[c][0] <= int(row['domain_a_end']) and int(row['domain_b_start']) <= mb[c][0] <= int(row['domain_b_end'])]
        pa, pb = np.array([ma[c][0] for c in columns]), np.array([mb[c][0] for c in columns])
        x = np.array([coords[links[ka]['model_path']][v] for v in pa]); y = np.array([coords[links[kb]['model_path']][v] for v in pb])
        xc, yc = x-x.mean(0), y-y.mean(0)
        rotation, _ = Rotation.align_vectors(yc, xc)
        residuals = np.sum((rotation.apply(xc)-yc)**2, axis=1)
        dx, dy = x[indices], y[indices]; dcx, dcy = dx-dx.mean(0), dy-dy.mean(0)
        drot, _ = Rotation.align_vectors(dcy, dcx)
        own = float(np.sqrt(np.mean(np.sum((drot.apply(dcx)-dcy)**2, axis=1))))
        distx, disty = pdist(dx), pdist(dy)
        ia, ib = np.triu_indices(len(indices), 1); ppa, ppb = pa[indices], pb[indices]
        local = ((distx <= 15) | (disty <= 15)) & (abs(ppa[ia]-ppa[ib]) >= 3) & (abs(ppb[ia]-ppb[ib]) >= 3)
        computed = {'whole_marker_rmsd_angstrom': float(np.sqrt(residuals.mean())),
                    'domain_sites_under_whole_marker_fit_rmsd_angstrom': float(np.sqrt(residuals[indices].mean())),
                    'domain_own_fit_rmsd_angstrom': own,
                    'local_mean_absolute_distance_change_angstrom': float(np.abs(distx-disty)[local].mean()) if local.any() else 0.0}
        if len(indices) != int(row['domain_compared_residues']) or int(local.sum()) != int(row['local_distance_pairs']):
            raise ValueError('Domain mask or local pair count mismatch')
        for key, value in computed.items():
            if not np.isclose(value, float(row[key]), atol=1e-8, rtol=1e-8):
                raise ValueError('Independent geometry mismatch: ' + key)
        records.append({k: row[k] for k in ['marker', 'taxon_a', 'taxon_b', 'pfam_accession']} | {'selection_sha256': digest, 'domain_residues': len(indices)} | computed)
    a.output.mkdir(parents=True); write_table(a.output / 'recalculated_domains.tsv', records)
    result = {'status': 'passed_independent_selected_domain_geometry', 'source_receipt_sha256': sha(a.domains / 'receipt.json'),
              'script_sha256': sha(Path(__file__)), 'comparisons': len(records),
              'selection': 'One smallest identity hash per Pfam accession, independent of metric values',
              'interpretation': 'Independent SciPy rotation and condensed-distance calculations reproduce selected whole/domain fits and local pair metrics. Only selected rows checked, not all domain geometry or PAE-filtered distances.',
              'artifacts': {'recalculated_domains.tsv': sha(a.output / 'recalculated_domains.tsv')}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n'); print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
