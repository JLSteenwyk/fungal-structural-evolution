#!/usr/bin/env python3
"""Superpose existing structures using fixed focal sequence-alignment mappings."""
import argparse
import hashlib
import json
from pathlib import Path
import sqlite3
import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio import SeqIO
from inspect_focal_domain_architectures import sha
from extract_domain_coordinates import load_atoms


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['mappings', 'sequences', 'bridge', 'output']:
        ap.add_argument('--' + name, type=Path, required=True)
    a = ap.parse_args()
    paths = [a.mappings, a.sequences, a.bridge, a.bridge.parent / 'receipt.json', Path(__file__)]
    paths += [Path(__file__).with_name(n + '.py') for n in
              ['inspect_focal_domain_architectures', 'extract_domain_coordinates', 'catalog_whole_proteome_structures']]
    pins = {str(p): sha(p) for p in paths}
    assert json.loads(paths[3].read_text())['artifacts'][a.bridge.name] == pins[str(a.bridge)]
    mappings = json.loads(a.mappings.read_text())
    assert mappings['status'] == 'complete_focal_alignment_sensitivity_mapping'
    assert mappings['source_hashes'][str(a.sequences)] == pins[str(a.sequences)]
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(a.sequences, 'fasta')}
    con = sqlite3.connect(a.bridge.resolve().as_uri() + '?mode=ro', uri=True)
    con.row_factory = sqlite3.Row
    models = {}; coordinates = {}; confidences = {}; missing = []
    for label, sequence in seqs.items():
        taxon, protein = label.split('_', 1)
        rows = con.execute('SELECT * FROM structures WHERE taxon_id=? AND protein_id=?', (taxon, protein)).fetchall()
        assert len(rows) <= 1
        if not rows:
            missing.append(label)
            continue
        row = dict(rows[0]); digest = hashlib.sha256(sequence.encode()).hexdigest()
        assert digest == row['sequence_sha256']
        path = row['model_path']; pins[path] = sha(path)
        recovered, atoms, ca = load_atoms(dict(path=path, sha256=pins[path], length=len(sequence), sequence_sha256=digest))
        assert recovered == sequence
        coordinates[label] = {p: np.array(next(v[2:5] for v in values if v[0] == 'CA')) for p, values in atoms.items()}
        confidences[label] = ca; models[label] = row
    con.close()
    results = []; skipped = []
    for mapping in mappings['results']:
        left, right = mapping['focal'], mapping['sister']
        key = {k: mapping[k] for k in ['mode', 'focal', 'sister', 'boundary', 'hit_id', 'domain_residues']}
        if left not in models or right not in models:
            skipped.append(key)
            continue
        for threshold in [0, 70, 90]:
            pairs = [(p, q) for p, q in mapping['residue_pairs']
                     if min(confidences[left][p], confidences[right][q]) >= threshold]
            rmsd = None
            if len(pairs) >= 3:
                x = np.array([coordinates[left][p] for p, q in pairs])
                y = np.array([coordinates[right][q] for p, q in pairs])
                fit = SVDSuperimposer(); fit.set(x, y); fit.run(); rmsd = float(fit.get_rms())
                assert np.linalg.det(fit.get_rotran()[0]) > 0
            results.append(dict(key, minimum_joint_plddt=threshold, matched_residues=len(pairs),
                                mapped_residues_before_confidence=len(mapping['residue_pairs']),
                                ca_rmsd_angstrom=rmsd, residue_pairs=pairs))
    assert all(sha(p) == h for p, h in pins.items())
    out = dict(status='complete_focal_fixed_mapping_geometry', source_hashes=pins, models=models,
               missing_from_frozen_bridge=missing, skipped=skipped, results=results,
               scope='Proper-rotation least-squares C-alpha RMSD using fixed sequence-alignment '
               'correspondences and joint pLDDT thresholds. No PAE filter, structural realignment, '
               'domain homology validation, or gain/loss inference. Missing bridge entries are not public-database absence.')
    with a.output.open('x') as f:
        json.dump(out, f, indent=2); f.write('\n')
    print(json.dumps([{k: v for k, v in r.items() if k != 'residue_pairs'} for r in results], indent=2))


if __name__ == '__main__':
    main()
