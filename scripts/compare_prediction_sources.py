#!/usr/bin/env python3
"""Compare every cross-pipeline exact-sequence model pair in a frozen inventory."""
import argparse
import csv
import fcntl
import gzip
import hashlib
import io
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from compare_marker_structures import ROOT, sha, geometry
from retrieve_matched_models import polymer_sequences
from retrieve_marker_pae import retrieve, validate_pae
from assess_pae_sensitivity import confidence_mask


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable output')
    raw = args.inventory.read_bytes()
    latest = {r['uniprot_accession']: r for r in map(json.loads, raw.decode().splitlines())}
    unique = {}
    for r in latest.values():
        if r['status'] == 'verified':
            for m in r['models']:
                unique[m['path']] = m
    groups = defaultdict(list)
    for m in unique.values():
        groups[m['sequence_sha256']].append(m)
    pairs = [(a, b) for models in groups.values() for a, b in combinations(sorted(models, key=lambda m: m['path']), 2)
             if (a.get('provider'), a.get('tool')) != (b.get('provider'), b.get('tool'))]
    args.output.mkdir(parents=True)
    config = {'inventory_path': str(args.inventory), 'inventory_sha256': hashlib.sha256(raw).hexdigest(),
        'distinct_models': len(unique), 'cross_pipeline_pairs': len(pairs),
        'script_sha256': sha(Path(__file__)),
        'resource_plan': 'One sequential CPU worker; reuse cached PAE, download missing PAE serially. Current frozen inventory has one pair of 173-residue models; less than 1 GB RAM and 10 MB output expected, seconds to minutes. No paid resources.'}
    (args.output / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
    cache = ROOT / 'data/structures/pae'
    lock = (cache / '.retrieval.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    loaded = {}
    provenance = {}
    for model in [m for pair in pairs for m in pair]:
        path = ROOT / model['path']
        if model['path'] in loaded:
            continue
        if sha(path) != model['sha256']:
            raise ValueError('Changed coordinates')
        cif = MMCIF2Dict(io.StringIO(path.read_text()))
        sequences = polymer_sequences(cif)
        if len(sequences) != 1 or hashlib.sha256(sequences[0].encode()).hexdigest() != model['sequence_sha256']:
            raise ValueError('Model polymer sequence differs')
        ca = {}
        for atom, p, x, y, z, confidence in zip(cif['_atom_site.label_atom_id'], cif['_atom_site.label_seq_id'],
                cif['_atom_site.Cartn_x'], cif['_atom_site.Cartn_y'], cif['_atom_site.Cartn_z'], cif['_atom_site.B_iso_or_equiv']):
            if atom == 'CA':
                if int(p) in ca:
                    raise ValueError('Ambiguous CA positions')
                ca[int(p)] = [float(x), float(y), float(z), float(confidence)]
        if set(ca) != set(range(1, len(sequences[0]) + 1)):
            raise ValueError('Incomplete CA positions')
        array = np.array([ca[p] for p in sorted(ca)])
        if not np.isfinite(array).all() or (array[:, 3] < 0).any() or (array[:, 3] > 100).any():
            raise ValueError('Invalid coordinates or confidence')
        pae_receipt = retrieve(model, cache)
        pae = validate_pae(gzip.decompress((ROOT / pae_receipt['path']).read_bytes()), model['length'])
        loaded[model['path']] = array, pae
        provenance[model['path']] = {'model': model, 'pae': pae_receipt}
    rows = []
    for a, b in pairs:
        aa, pae_a = loaded[a['path']]; bb, pae_b = loaded[b['path']]
        if aa.shape != bb.shape:
            raise ValueError('Exact-sequence model lengths differ')
        for cutoff in [0, 70, 90]:
            mask = (aa[:, 3] >= cutoff) & (bb[:, 3] >= cutoff)
            positions = np.flatnonzero(mask) + 1
            row = {'sequence_sha256': a['sequence_sha256'], 'model_a': a['model_id'], 'version_a': a['version'],
                'provider_a': a.get('provider'), 'tool_a': a.get('tool'), 'model_b': b['model_id'], 'version_b': b['version'],
                'provider_b': b.get('provider'), 'tool_b': b.get('tool'), 'protein_length': len(aa),
                'plddt_cutoff': cutoff, 'matched_residues': len(positions), 'retained_fraction': len(positions) / len(aa),
                'status': 'compared' if len(positions) >= 50 and len(positions) >= .5 * len(aa) else 'insufficient_coverage'}
            metrics = dict.fromkeys(['ca_superposition_rmsd_angstrom', 'local_distance_pairs', 'local_distance_mean_absolute_change_angstrom',
                'local_distance_rms_change_angstrom', 'pae10_local_pairs', 'pae10_local_retained_fraction', 'pae10_local_mean_absolute_change_angstrom'], '')
            if row['status'] == 'compared':
                x, y = aa[mask, :3], bb[mask, :3]
                metrics.update(geometry(x, y, positions, positions))
                dx = np.linalg.norm(x[:, None] - x[None, :], axis=2)
                dy = np.linalg.norm(y[:, None] - y[None, :], axis=2)
                local = (np.triu(np.ones(dx.shape, dtype=bool), 1) & ((dx <= 15) | (dy <= 15))
                         & (np.abs(positions[:, None] - positions[None, :]) >= 3))
                retained = local & confidence_mask(pae_a, pae_b, positions, positions, 10)
                metrics.update(pae10_local_pairs=int(retained.sum()),
                    pae10_local_retained_fraction=float(retained.sum() / local.sum()) if local.any() else '',
                    pae10_local_mean_absolute_change_angstrom=float(np.abs(dx[retained] - dy[retained]).mean()) if retained.any() else '')
            rows.append(row | metrics)
    if rows:
        with (args.output / 'source_comparisons.tsv').open('w') as handle:
            w = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
    (args.output / 'model_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    receipt = {'status': 'complete_available_cross_pipeline_comparison', 'config_sha256': sha(args.output / 'config.json'),
        'pairs': len(pairs), 'unique_sequences': len({a['sequence_sha256'] for a, b in pairs}), 'comparison_rows': len(rows),
        'interpretation': 'Every cross-pipeline exact-sequence pair available in the frozen inventory. Confidence cutoffs condition the comparison and do not calibrate accuracy across predictors. A single available sequence is an initial control within the full atlas, not evidence for general systematic source effects. Both pipelines derive structures from sequence; experimental validation and broader alternative-method comparisons remain necessary.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    main()
