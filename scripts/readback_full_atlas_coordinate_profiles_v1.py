#!/usr/bin/env python3
"""Independently decode every original source and every exported full-atlas C-alpha profile."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import multiprocessing
from pathlib import Path
import tarfile

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.PDBData import protein_letters_3to1

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def independent_arrays(model):
    data = Path(model['path']).read_bytes()
    if hashlib.sha256(data).hexdigest() != model['sha256']:
        raise RuntimeError('Immutable source coordinate bytes changed')
    fields = MMCIF2Dict(io.StringIO(data.decode()))
    polymers = [''.join(s.split()) for s in fields['_entity_poly.pdbx_seq_one_letter_code_can']]
    if len(polymers) != 1:
        raise ValueError('Multiple polymer entities')
    sequence = polymers[0]
    if len(sequence) != model['length'] or hashlib.sha256(sequence.encode()).hexdigest() != model['sequence_sha256']:
        raise ValueError('Source sequence identity differs')
    keys = ['group_PDB', 'label_atom_id', 'label_comp_id', 'label_seq_id', 'label_asym_id',
            'label_alt_id', 'pdbx_PDB_model_num', 'Cartn_x', 'Cartn_y', 'Cartn_z',
            'occupancy', 'B_iso_or_equiv', 'type_symbol']
    columns = {key: fields['_atom_site.'+key] for key in keys}
    if len({len(v) for v in columns.values()}) != 1:
        raise ValueError('Atom columns differ in length')
    seen, chains, atoms, ca = set(), set(), {}, {}
    for index in range(len(columns['label_atom_id'])):
        value = {key: columns[key][index] for key in keys}
        if value['group_PDB'] != 'ATOM' or value['label_alt_id'] not in ('.', '?') or value['pdbx_PDB_model_num'] != '1':
            raise ValueError('Unsupported source atom/model/alternate location')
        position, name = int(value['label_seq_id']), value['label_atom_id']
        if not 1 <= position <= len(sequence) or protein_letters_3to1.get(value['label_comp_id']) != sequence[position-1]:
            raise ValueError('Noncanonical or differing residue identity')
        if (position, name) in seen:
            raise ValueError('Duplicate atom')
        seen.add((position, name)); chains.add(value['label_asym_id'])
        numbers = [float(value[key]) for key in ('Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv')]
        if not all(math.isfinite(v) for v in numbers) or not 0 <= numbers[3] <= 1 or not 0 <= numbers[4] <= 100:
            raise ValueError('Invalid source coordinate, occupancy or confidence')
        if len(name) > 4 or len(value['type_symbol']) > 2:
            raise ValueError('Unsupported atom/element label')
        atoms.setdefault(position, set()).add(name)
        if name == 'CA':
            ca[position] = numbers[:3]+[numbers[4]]
    if len(chains) != 1 or sorted(ca) != list(range(1, len(sequence)+1)):
        raise ValueError('C-alpha coverage/chain identity differs')
    matrix = np.asarray([ca[p] for p in range(1, len(sequence)+1)], dtype=np.float64)
    return dict(sequence_ascii=np.frombuffer(sequence.encode('ascii'), dtype=np.uint8).copy(),
                ca_xyz=matrix[:, :3].copy(), ca_plddt=matrix[:, 3].copy(),
                missing_backbone=np.asarray([bool({'N', 'CA', 'C', 'O'}-atoms[p])
                                            for p in range(1, len(sequence)+1)], dtype=np.bool_))


def check_profile(blob, row, expected, model):
    assert hashlib.sha256(blob).hexdigest() == row['npz_sha256']
    with np.load(io.BytesIO(blob), allow_pickle=False) as actual:
        assert set(actual.files) == set(expected)
        for key, value in expected.items():
            assert actual[key].dtype == value.dtype and actual[key].shape == value.shape
            assert np.array_equal(actual[key], value), 'Original source array differs: '+key
            digest = hashlib.sha256(value.dtype.str.encode()+str(value.shape).encode()+value.tobytes()).hexdigest()
            assert row['array_sha256'][key] == digest
    confidence = expected['ca_plddt']
    n = len(confidence)
    for cutoff in (50, 70, 90):
        # Scalar accumulation is independent of the producer's vectorized masks.
        assert row['plddt_ge'+str(cutoff)] == sum(float(v) >= cutoff for v in confidence)
    assert row['residues'] == n == model['length']
    assert row['missing_backbone_residues'] == sum(bool(v) for v in expected['missing_backbone'])
    assert row['minimum_ca_plddt'] == min(map(float, confidence))
    assert row['maximum_ca_plddt'] == max(map(float, confidence))
    assert row['catalog_mean_ca_plddt'] == model['mean_ca_plddt']
    mean = math.fsum(map(float, confidence))/n
    assert abs(row['mean_coordinate_ca_plddt']-mean) <= 1e-10
    assert abs(row['coordinate_minus_catalog_mean_ca_plddt']-(mean-model['mean_ca_plddt'])) <= 1e-10


def read_shard(proof, shard_dir):
    job = Path(proof['job']); directory = Path(shard_dir); r = proof['receipt']
    assert sha(job) == r['job_sha256']
    for name, digest in r['artifacts'].items():
        assert sha(directory/name) == digest
    counts = {s: Counter(dict.fromkeys(r['counts'][s], 0)) for s in ('AFDB', 'ESMFold')}
    with tarfile.open(directory/(job.stem+'.tar')) as archive, job.open() as inp, (directory/(job.stem+'.jsonl')).open() as rows:
        members = archive.getmembers()
        names = [m.name for m in members]
        assert len(set(names)) == len(names) and all(m.isfile() for m in members)
        remaining = set(names)
        emitted = (json.loads(line) for line in rows)
        for line in inp:
            raw = json.loads(line); source, model = raw['source'], raw['model']
            row = next(emitted)
            identity = dict(source=source, sequence_sha256=model['sequence_sha256'], model_id=model['model_id'],
                            version=model['version'], source_path=model['path'], source_sha256=model['sha256'],
                            length=model['length'], prediction_config_sha256=model.get('prediction_config_sha256', ''))
            assert all(row[key] == value for key, value in identity.items()) and row['pae_evaluated'] is False
            c = counts[source]; c['models'] += 1; c['source_coordinate_bytes'] += Path(model['path']).stat().st_size
            try:
                expected = independent_arrays(model)
            except (ValueError, KeyError):
                assert row['status'] == 'source_validation_rejected' and row.get('reason')
                assert 'member' not in row
                c['rejected_models'] += 1
            else:
                assert row['status'] == 'coordinate_and_residue_profile_exported'
                assert row['member'] == source+'-'+model['sequence_sha256']+'.npz'
                assert row['member'] in remaining
                remaining.remove(row['member'])
                check_profile(archive.extractfile(row['member']).read(), row, expected, model)
                c['valid_models'] += 1; c['valid_residues'] += row['residues']
                for key in ('missing_backbone_residues', 'plddt_ge50', 'plddt_ge70', 'plddt_ge90'):
                    c[key] += row[key]
        assert next(emitted, None) is None and not remaining
    assert {s: dict(c) for s, c in counts.items()} == r['counts']
    return {s: dict(c) for s, c in counts.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text()); pins = dict(plan['pins']); verify(pins)
    pp = Path(plan['producer_plan']); producer_plan = json.loads(pp.read_text())
    own = Path(plan['producer_receipt']); producer = json.loads(own.read_text())
    tp = Path(plan['producer_transport']); transport = json.loads(tp.read_text())
    assert producer['status'] == 'completed_full_atlas_coordinate_and_residue_audit_pending_independent_readback'
    assert transport['status'] == 'verified_original_software_wait_and_whole_wrapper_payloads'
    assert transport['validation_sha256'] == sha(own) and transport['original_tool_terminal_exit_code'] == 0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    assert transport['manager_start_records'] == transport['manager_completion_records'] == 1
    for mapping in (producer['source_hashes'], transport['source_hashes'], producer_plan['pins']):
        for path, digest in mapping.items():
            bind(pins, path, digest)
    verify(pins)
    root = Path(producer_plan['output']); raw_path = root/'receipt.json'; raw = json.loads(raw_path.read_text())
    assert sha(raw_path) == producer['raw_receipt_sha256'] and raw['counts'] == producer['counts']
    assert raw['plan_sha256'] == sha(pp)
    output = Path(plan['output']); assert not output.exists(); output.mkdir(parents=True)
    sources = Counter(); proofs = sorted(raw['proofs'], key=lambda p: p['job'])
    with Path(producer_plan['models']).open() as original:
        expected = (json.loads(line) for line in original)
        for proof in proofs:
            with Path(proof['job']).open() as jobs:
                for line in jobs:
                    observed = json.loads(line)
                    assert observed == next(expected), 'Full original source/job scope differs'
                    sources[observed['source']] += 1
        assert next(expected, None) is None
    assert dict(sources) == producer_plan['expected_models_by_source']
    counts = {s: Counter(dict.fromkeys(raw['counts'][s], 0)) for s in sources}; completed = 0
    with ProcessPoolExecutor(max_workers=plan['resources']['cpu'], mp_context=multiprocessing.get_context('spawn')) as pool:
        futures = [pool.submit(read_shard, proof, str(root/'shards')) for proof in proofs]
        for future in as_completed(futures):
            try:
                checked = future.result()
            except BaseException:
                for pending in futures:
                    pending.cancel()
                raise
            for source, local in checked.items():
                counts[source].update(local)
            completed += 1
            state = dict(stage='independent_readback_running', completed_shards=completed,
                         counts={s: dict(c) for s, c in counts.items()})
            temporary = output/'state.tmp'; temporary.write_text(json.dumps(state, indent=2)+'\n')
            temporary.replace(output/'state.json')
    assert {s: dict(c) for s, c in counts.items()} == raw['counts']
    assert completed == producer['shards'] == raw['shards']
    for path in (args.plan, pp, own, tp, raw_path, Path(__file__)):
        bind(pins, path)
    verify(pins)
    result = dict(status='passed_full_atlas_original_coordinate_profile_and_rejection_readback',
                  checked_utc=datetime.now(timezone.utc).isoformat(), counts=raw['counts'], shards=completed,
                  source_hashes=pins, scientific_eligibility=False, gpu=False, new_predictions=0,
                  scope='Every original model/job row, all source byte hashes/atom/sequence identities, '
                        'every profile value/dtype/shape/archive member and scalar threshold count '
                        'independently reconstructed without producer imports. Rejection eligibility '
                        'rederived; recorded exact error causation is not independently adjudicated. '
                        'Shared CIF lexical parser and canonical residue map; no PAE/context features, '
                        'prediction calibration, homology or evolutionary-event acceptance.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
