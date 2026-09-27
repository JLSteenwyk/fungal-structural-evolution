"""Compare identical coordinate subsets around all shared functional sites."""
import argparse
from collections import Counter
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.PDBData import protein_letters_3to1
from assess_reference_alignment_geometry import geometry
from readback_reference_alignment_geometry_v2 import verify_row, quaternion_curvature


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def coordinates(model):
    assert sha(model['path']) == model['sha256']
    raw = MMCIF2Dict(model['path'])
    fields = ['label_atom_id', 'label_seq_id', 'label_comp_id', 'Cartn_x', 'Cartn_y',
              'Cartn_z', 'B_iso_or_equiv', 'label_asym_id', 'pdbx_PDB_model_num']
    columns = [raw['_atom_site.' + key] for key in fields]
    assert len({len(c) for c in columns}) == 1
    residues, chains, models = {}, set(), set()
    for atom, pos, aa, x, y, z, confidence, chain, number in zip(*columns):
        chains.add(chain)
        models.add(number)
        if atom != 'CA':
            continue
        index = int(pos)
        assert index not in residues
        residues[index] = (protein_letters_3to1[aa], [float(x), float(y), float(z)], float(confidence))
    assert len(chains) == len(models) == 1
    assert set(residues) == set(range(1, model['length'] + 1))
    sequence = ''.join(residues[i][0] for i in sorted(residues))
    assert hashlib.sha256(sequence.encode()).hexdigest() == model['sequence_sha256']
    xyz = np.array([residues[i][1] for i in sorted(residues)])
    confidence = np.array([residues[i][2] for i in sorted(residues)])
    assert np.isfinite(xyz).all() and np.isfinite(confidence).all()
    assert ((confidence >= 0) & (confidence <= 100)).all()
    return sequence, xyz, confidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    source = Path('results/functional_sites/prediction-context-comparison-20260927-v1')
    receipt = json.loads((source / 'receipt.json').read_text())
    assert receipt['status'] == 'complete_full_joint_functional_context_comparison'
    table = source / 'functional_context_comparison.tsv'
    assert sha(table) == receipt['artifacts'][table.name]
    rows = list(csv.DictReader(table.open(), delimiter='\t'))
    assert len(rows) == receipt['residues'] == 150
    roots = dict(afdb=Path('results/structural_markers/afdb-recovered-union-20260925-v1/mapping'),
                 esmfold=Path('results/structural_markers/esmfold-all-completed-20260922-v1'))
    joins = dict(afdb=Path('results/functional_sites/afdb-recovered-linked-20260927-v1'),
                 esmfold=Path('results/functional_sites/esmfold-all-completed-linked-20260923-v1'))
    pins = {str(source / 'receipt.json'): sha(source / 'receipt.json'), str(table): sha(table)}
    provenance, cache = {}, {}
    for label, root in roots.items():
        join = json.loads((joins[label] / 'receipt.json').read_text())
        assert receipt['source_sha256'][str(joins[label] / 'receipt.json')] == sha(joins[label] / 'receipt.json')
        assert join['source_receipts']['snapshot']['sha256'] == sha(root / 'receipt.json')
        rp = json.loads((root / 'receipt.json').read_text())
        mp = root / 'model_provenance.json'
        assert sha(mp) == rp['artifacts'][mp.name]
        entries = json.loads(mp.read_text())
        provenance[label] = {(r['model_id'], str(r['version'])): r for r in entries}
        assert len(entries) == len(provenance[label])
        for path in (mp, root / 'receipt.json', joins[label] / 'receipt.json'):
            pins[str(path)] = sha(path)
    output, matrices = [], {}
    maximum_squared_error = 0.
    for r in rows:
        models = []
        for label in roots:
            key = label, r[label + '_model_id'], r[label + '_model_version']
            model = provenance[label][key[1:]]
            assert model['sequence_sha256'] == r['sequence_sha256']
            if key not in cache:
                cache[key] = coordinates(model)
                pins[model['path']] = model['sha256']
            models.append(cache[key])
            own_context = np.array(list(map(int, r[label + '_context_residues_1based'].split(',')))) - 1
            assert math.isclose(float(min(cache[key][2][own_context])), float(r[label + '_feature_min_plddt']), abs_tol=1e-10)
        assert models[0][0] == models[1][0]
        pos = int(r['protein_residue_1based'])
        a = set(map(int, r['afdb_context_residues_1based'].split(',')))
        e = set(map(int, r['esmfold_context_residues_1based'].split(',')))
        definitions = dict(focal_triplet={pos-1, pos, pos+1}, afdb_context=a, esmfold_context=e, context_union=a | e)
        for definition, positions in definitions.items():
            indices = np.array(sorted(positions)) - 1
            x, y = [m[1][indices] for m in models]
            g = geometry(x, y)
            error, near_zero = verify_row(g, x, y)
            xc, yc = x - x.mean(0), y - y.mean(0)
            u, _, vt = np.linalg.svd(xc.T @ yc)
            correction = np.diag([1., 1., 1. if np.linalg.det(u @ vt) >= 0 else -1.])
            rmsd = float(np.sqrt(np.mean(np.sum((xc @ u @ correction @ vt - yc)**2, axis=1))))
            _, optimum = quaternion_curvature(x, y)
            squared = float((np.sum(xc**2) + np.sum(yc**2) - 2*optimum) / len(x))
            assert math.isclose(rmsd**2, squared, rel_tol=1e-9, abs_tol=1e-9)
            maximum_squared_error = max(maximum_squared_error, abs(rmsd**2 - squared))
            differences = [math.dist(x[i], x[j]) - math.dist(y[i], y[j]) for i, j in itertools.combinations(range(len(x)), 2)]
            d = np.linalg.norm(x[:, None] - x[None, :], axis=2) - np.linalg.norm(y[:, None] - y[None, :], axis=2)
            np.testing.assert_allclose(d[np.triu_indices(len(x), 1)], differences, atol=1e-12, rtol=1e-12)
            record = {k: r[k] for k in ['marker', 'taxon_id', 'protein_id', 'protein_residue_1based',
                'paired_state_comparison', 'partner_comparison', 'any_conserved_candidate']}
            record.update(context_definition=definition, residue_positions_1based=','.join(map(str, sorted(positions))),
                compared_residues=len(x), ca_rmsd_angstrom=rmsd, quaternion_rmsd_squared=squared,
                pair_distance_comparisons=len(differences), pair_distance_mean_absolute_difference_angstrom=float(np.mean(np.abs(differences))),
                pair_distance_rms_difference_angstrom=float(np.sqrt(np.mean(np.square(differences)))),
                minimum_ca_plddt_afdb=float(min(models[0][2][indices])), minimum_ca_plddt_esmfold=float(min(models[1][2][indices])),
                rank_afdb=g['rank_left'], rank_esmfold=g['rank_right'], rotation_status=g['geometry_status'],
                quaternion_scaled_curvature_error=error, near_zero_quaternion_gap=near_zero)
            output.append(record)
            matrices[tuple(record[k] for k in ['marker', 'taxon_id', 'protein_residue_1based', 'context_definition'])] = (x, y)
    assert len(output) == 600 and len(matrices) == 600
    args.output.mkdir(parents=True)
    path = args.output / 'functional_context_geometry.tsv'
    with path.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(output)
    exported = list(csv.DictReader(path.open(), delimiter='\t'))
    for r in exported:
        x, y = matrices[tuple(r[k] for k in ['marker', 'taxon_id', 'protein_residue_1based', 'context_definition'])]
        _, optimum = quaternion_curvature(x, y)
        squared = (np.sum((x-x.mean(0))**2) + np.sum((y-y.mean(0))**2) - 2*optimum) / len(x)
        assert math.isclose(float(r['ca_rmsd_angstrom'])**2, squared, rel_tol=1e-9, abs_tol=1e-9)
        distances = [abs(math.dist(x[i], x[j])-math.dist(y[i], y[j])) for i, j in itertools.combinations(range(len(x)), 2)]
        assert math.isclose(float(r['pair_distance_mean_absolute_difference_angstrom']), sum(distances)/len(distances), abs_tol=1e-12)
    for p, h in pins.items():
        assert sha(p) == h
    result = dict(status='complete_functional_context_coordinate_comparison', focal_residues=150, comparison_rows=600,
        coordinate_models=len(cache), rotation_status_counts=dict(Counter(r['rotation_status'] for r in output)),
        maximum_svd_quaternion_squared_rmsd_error=maximum_squared_error,
        source_sha256=pins, script_sha256=sha(__file__), artifacts={path.name: sha(path)},
        scope='Four exact sequence-position subsets for every shared functional residue; raw single-chain CA sequences and hashes checked. Unique positions receive equal weight. SVD RMSD/curvature cross-checked with quaternion optimum; scalar and vector pair distances agree. Cross-predictor contexts retain confidence without additional filtering. Serialized RMSD and mean-distance values read back. Not whole-protein displacement, native descriptor reconstruction, predictor accuracy, side-chain geometry, functional validation or evolutionary change.')
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256','artifacts')}, indent=2))


if __name__ == '__main__':
    main()
