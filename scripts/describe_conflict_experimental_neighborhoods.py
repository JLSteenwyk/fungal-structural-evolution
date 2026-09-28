#!/usr/bin/env python3
"""Describe deposited heavy-atom neighborhoods of experimentally mapped conflict sites."""
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    base = Path('results/ancestral/ancestral-conflict-experimental-coverage-20260927-v1')
    source = base / 'conflict_experimental_coverage.json'
    assert sha(source) == json.loads((base / 'receipt.json').read_text())['artifacts'][source.name]
    pins = {str(p): sha(p) for p in [source, base / 'receipt.json']}
    targets = {}
    for row in json.loads(source.read_text()):
        for context in row['experimental_contexts']:
            for observed in context['observed_positions']:
                key = tuple(observed[k] for k in ['entry_id', 'model_number', 'label_asym_id', 'label_seq_id'])
                if key in targets:
                    assert targets[key] == observed
                targets[key] = observed
    assert len(targets) == 8
    coordinates = Path('results/experimental_structures/whole-domain-case-coordinates-20260927-v1')
    crp = coordinates / 'receipt.json'
    archives = {r['entry_id']: r for r in json.loads(crp.read_text())['results']}
    config = Path('results/experimental_structures/whole-domain-case-ca-mapping-20260927-v1/config.json')
    assert sha(crp) == json.loads(config.read_text())['coordinate_receipt_sha256']
    pins.update({str(p): sha(p) for p in [crp, config]})
    rows = []
    for entry in sorted({k[0] for k in targets}):
        path = coordinates / (entry + '.cif.gz')
        assert sha(path) == archives[entry]['gzip_sha256']
        pins[str(path)] = sha(path)
        with gzip.open(path, 'rt') as stream:
            data = MMCIF2Dict(stream)
        types = dict(zip(data['_entity.id'], data['_entity.type']))
        names = ['id', 'type_symbol', 'label_atom_id', 'label_entity_id', 'label_asym_id',
                 'label_seq_id', 'label_comp_id', 'label_alt_id', 'pdbx_PDB_model_num',
                 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code']
        columns = [data['_atom_site.' + name] for name in names]
        assert len({len(c) for c in columns}) == 1
        atoms = [dict(zip(names, values)) for values in zip(*columns)]
        atoms = [a for a in atoms if a['type_symbol'].upper() not in ['H', 'D'] and float(a['occupancy']) > 0]
        for a in atoms:
            a['xyz'] = [float(a['Cartn_' + axis]) for axis in 'xyz']
            assert all(math.isfinite(x) for x in a['xyz'])
        for key, observed in targets.items():
            if key[0] != entry:
                continue
            _, model, chain, position = key
            focal = [a for a in atoms if (a['pdbx_PDB_model_num'], a['label_asym_id'], a['label_seq_id']) == (model, chain, position)]
            assert focal and all(a['label_comp_id'] == 'TYR' for a in focal)
            ca = [a for a in focal if a['label_atom_id'] == 'CA']
            assert len(ca) == 1
            prior_ca = json.loads(observed['atom_records_json'])
            assert len(prior_ca) == 1
            assert ca[0]['xyz'] == [float(prior_ca[0]['Cartn_' + axis]) for axis in 'xyz']
            # Retain all positive-occupancy alternate records; no conformation is selected.
            residues = {}
            for atom in atoms:
                if atom['pdbx_PDB_model_num'] != model or atom in focal:
                    continue
                etype = types[atom['label_entity_id']]
                if etype == 'polymer':
                    if atom['label_asym_id'] == chain:
                        if abs(int(atom['label_seq_id']) - int(position)) <= 1:
                            continue
                        category = 'same_chain_nonadjacent_polymer'
                    else:
                        category = 'other_chain_polymer'
                elif etype == 'water':
                    category = 'water'
                else:
                    category = 'nonpolymer_or_branched'
                distances = np.linalg.norm(np.asarray([a['xyz'] for a in focal]) - atom['xyz'], axis=1)
                i = int(np.argmin(distances))
                distance = float(distances[i])
                assert math.isclose(distance, min(math.dist(a['xyz'], atom['xyz']) for a in focal), abs_tol=1e-10)
                identity = (category, atom['label_asym_id'], atom['label_seq_id'], atom['auth_seq_id'],
                            atom['pdbx_PDB_ins_code'], atom['label_comp_id'])
                if identity not in residues or distance < residues[identity]['distance_angstrom']:
                    residues[identity] = dict(category=category, chain=identity[1], label_seq_id=identity[2],
                        auth_seq_id=identity[3], insertion_code=identity[4], component=identity[5],
                        distance_angstrom=distance, focal_atom=focal[i], neighbor_atom=atom)
            nearest = {category: min((r for r in residues.values() if r['category'] == category),
                                     key=lambda r: r['distance_angstrom'], default=None)
                       for category in ['same_chain_nonadjacent_polymer', 'other_chain_polymer', 'water', 'nonpolymer_or_branched']}
            contacts = sorted((r for r in residues.values() if r['distance_angstrom'] <= 5), key=lambda r:r['distance_angstrom'])
            rows.append(dict(entry=entry, model=model, chain=chain, entity_position=int(position),
                residue='Y', focal_atom_names=sorted({a['label_atom_id'] for a in focal}),
                focal_records=focal, nearest_by_category=nearest, residue_contacts_within_5_angstrom=contacts))
    out = Path('results/ancestral/ancestral-conflict-experimental-neighborhoods-20260927-v1')
    out.mkdir(exist_ok=False)
    path = out / 'deposited_neighborhoods.json'
    path.write_text(json.dumps(rows, indent=2) + '\n')
    result = dict(status='complete_deposited_conflict_site_neighborhoods', sites=len(rows), entries=2,
        pins=pins, script_sha256=sha(__file__), artifacts={path.name:sha(path)},
        verification='Every minimum atom distance checked using independent scalar math.dist; focal CA coordinates match the prior correspondence table exactly.',
        scope='Deposited coordinates only, no assembly or symmetry expansion. All positive-occupancy heavy-atom alternatives retained; minima may use different alternate conformers. Same-chain sequence-adjacent residues excluded. Contacts are geometric descriptions, not validated interactions, binding sites, catalytic roles, ancestral polarity or independent experimental replicates.')
    rp = out / 'receipt.json'
    rp.write_text(json.dumps(result, indent=2) + '\n')
    result.update(completed_receipt_path=str(rp), completed_receipt_sha256=sha(rp))
    Path('metadata/ancestral_conflict_experimental_neighborhoods_completed_20260927.json').write_text(json.dumps(result, indent=2) + '\n')
    for row in rows:
        print(row['entry'], row['chain'], {k: None if v is None else (v['component'], round(v['distance_angstrom'], 3)) for k,v in row['nearest_by_category'].items()})


if __name__ == '__main__':
    main()
