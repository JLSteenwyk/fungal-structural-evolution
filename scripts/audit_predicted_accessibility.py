#!/usr/bin/env python3
"""Audit completed ASA entries against raw atom tables; never infer biological exposure."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.IUPACData import protein_letters_3to1
from audit_busco_gene_copies import ROOT, sha
from assess_pae_sensitivity import checked_receipt
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['assessment', 'snapshot', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--allow-partial', action='store_true')
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable audit output')
    checked_receipt(a.snapshot)
    cp = a.assessment / 'config.json'
    c = json.loads(cp.read_text()); config_hash = sha(cp)
    if c['snapshot_receipt_sha256'] != sha(a.snapshot / 'receipt.json'):
        raise ValueError('Snapshot differs')
    models = json.loads((a.snapshot / 'model_provenance.json').read_text())
    ids = {m['model_id'] for m in models}
    if len(ids) != len(models) or len(models) != c['models'] or sum(m['length'] for m in models) != c['residues']:
        raise ValueError('Configured model grid differs')
    # Freeze the available atomic entry receipts before any expensive reads.
    entries = {f.name.removesuffix('.receipt.json'): sha(f) for f in a.assessment.glob('*.receipt.json')}
    if not set(entries) <= ids or (not a.allow_partial and set(entries) != ids):
        raise ValueError('Missing or unexpected completed entries')
    if not entries:
        raise ValueError('No completed entries')
    full_receipt = a.assessment / 'receipt.json'
    full_hash = sha(full_receipt) if full_receipt.exists() else None
    if not a.allow_partial and full_hash is None:
        raise ValueError('Full production receipt required')
    if full_hash:
        r = json.loads(full_receipt.read_text())
        if r['config_sha256'] != config_hash or r['entry_receipts'] != entries or r['models'] != len(models) or r['residues'] != c['residues']:
            raise ValueError('Full completion receipt differs')
    aa = {k.upper(): v for k, v in protein_letters_3to1.items()}
    summaries = []
    for m in models:
        sid = m['model_id']
        if sid not in entries:
            continue
        rp = a.assessment / (sid + '.receipt.json')
        r = json.loads(rp.read_text()); table = a.assessment / (sid + '.residues.tsv.gz')
        if sha(rp) != entries[sid] or r['model_id'] != sid or r['config_sha256'] != config_hash or r['model_sha256'] != m['sha256'] or r['sequence_sha256'] != m['sequence_sha256'] or sha(table) != r['table_sha256']:
            raise ValueError('Entry provenance differs: ' + sid)
        path = ROOT / m['path']
        if sha(path) != m['sha256']:
            raise ValueError('Coordinates changed')
        raw = MMCIF2Dict(str(path)); atoms = defaultdict(list)
        fields = ['label_seq_id', 'label_comp_id', 'label_atom_id', 'type_symbol', 'B_iso_or_equiv', 'label_asym_id', 'pdbx_PDB_model_num']
        columns = [raw['_atom_site.' + k] for k in fields]
        if len({len(x) for x in columns}) != 1:
            raise ValueError('Atom columns differ')
        chains = set(); model_numbers = set()
        for pos, res, atom, elem, b, chain, number in zip(*columns):
            if elem not in {'C', 'N', 'O', 'S'}:
                raise ValueError('Unexpected atom element')
            atoms[int(pos)].append((res, atom, float(b)))
            chains.add(chain); model_numbers.add(number)
        if len(chains) != 1 or len(model_numbers) != 1 or set(atoms) != set(range(1, m['length'] + 1)):
            raise ValueError('Raw coordinate grid differs')
        with gzip.open(table, 'rt') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        if len(rows) != m['length'] or len(rows) != r['residues']:
            raise ValueError('Residue count differs')
        areas = []; sequence = []; total_atoms = 0
        for i, row in enumerate(rows, 1):
            at = atoms[i]; residues = {x[0] for x in at}; ca = [x for x in at if x[1] == 'CA']
            if len(residues) != 1 or len(ca) != 1 or len({x[1] for x in at}) != len(at):
                raise ValueError('Ambiguous raw residue')
            residue_aa = aa[next(iter(residues))]; sequence.append(residue_aa)
            if row['model_id'] != sid or row['sequence_sha256'] != m['sequence_sha256'] or int(row['protein_residue_1based']) != i or row['amino_acid'] != residue_aa or int(row['heavy_atom_count']) != len(at) or row['context'] != 'isolated_predicted_chain':
                raise ValueError('Residue identity/atom count differs')
            confidence = float(row['ca_plddt']); area = float(row['sasa_angstrom_squared'])
            if not math.isfinite(confidence) or not 0 <= confidence <= 100 or confidence != ca[0][2] or not math.isfinite(area) or area < 0:
                raise ValueError('Invalid ASA/confidence')
            areas.append(area); total_atoms += len(at)
        if hashlib.sha256(''.join(sequence).encode()).hexdigest() != m['sequence_sha256'] or total_atoms != r['heavy_atoms'] or not math.isclose(math.fsum(areas), r['total_sasa_angstrom_squared'], rel_tol=1e-12, abs_tol=1e-8):
            raise ValueError('Entry totals differ')
        summaries.append({'model_id': sid, 'residues': len(rows), 'heavy_atoms': total_atoms, 'total_sasa_angstrom_squared': math.fsum(areas), 'entry_receipt_sha256': entries[sid]})
        if len(summaries) % 250 == 0:
            print('audited', len(summaries), '/', len(entries), flush=True)
    if sha(cp) != config_hash:
        raise ValueError('Configuration changed during audit')
    a.output.mkdir(parents=True)
    write_table(a.output / 'audited_models.tsv', summaries)
    out = {'status': 'passed_frozen_completed_entries' if full_hash is None else 'passed_full_accessibility_snapshot',
           'snapshot_receipt_sha256': sha(a.snapshot / 'receipt.json'), 'assessment_config_sha256': config_hash,
           'assessment_receipt_sha256': full_hash, 'script_sha256': sha(Path(__file__)),
           'models_audited': len(summaries), 'models_expected': len(models), 'residues_audited': sum(x['residues'] for x in summaries),
           'scope': 'All frozen completed entries checked against raw mmCIF atom tables, including hashes, sequence, residue grid, atom counts, CA confidence and ASA totals. ASA itself was not independently recalculated. Partial completion does not validate unaudited or unfinished entries. No biological exposure, core, surface or interface validation.',
           'artifacts': {'audited_models.tsv': sha(a.output / 'audited_models.tsv')}}
    (a.output / 'receipt.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
