#!/usr/bin/env python3
"""Run separately tracked ESMFold controls with configurable attention chunking."""
import argparse
import fcntl
import hashlib
import importlib.metadata
import json
import os
import platform
import signal
import time
from pathlib import Path
import numpy as np
import torch
from transformers import EsmForProteinFolding
from transformers.models.esm import modeling_esmfold
from prepare_pfam import ROOT, digest

AA3 = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),
               'ARNDCQEGHILKMFPSTWYV'))


def read_fasta(path):
    sid, parts = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if sid is not None:
                    yield sid, ''.join(parts)
                sid, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
    if sid is not None:
        yield sid, ''.join(parts)


def validate_prediction(sequence, pdb, pae, ca_plddt, max_pae):
    ca = [line for line in pdb.splitlines() if line.startswith('ATOM') and line[12:16].strip() == 'CA']
    if (''.join(AA3.get(line[17:20], '?') for line in ca) != sequence
            or [int(line[22:26]) for line in ca] != list(range(1, len(sequence) + 1))
            or len({line[21] for line in ca}) != 1):
        raise ValueError('Predicted PDB does not preserve complete sequence/residue numbering')
    xyz = np.array([[float(line[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)]] for line in ca])
    pdb_confidence = np.array([float(line[60:66]) for line in ca])
    if (not np.isfinite(xyz).all() or ca_plddt.shape != (len(sequence),)
            or not np.isfinite(ca_plddt).all() or np.any(ca_plddt < 0) or np.any(ca_plddt > 100)
            or not np.allclose(pdb_confidence, ca_plddt, atol=.0051, rtol=0)):
        raise ValueError('Invalid coordinates or incorrectly scaled PDB confidence')
    if (pae.shape != (len(sequence), len(sequence)) or not np.isfinite(pae).all()
            or np.any(pae < 0) or not np.isfinite(max_pae) or max_pae <= 0
            or np.any(pae > max_pae + 1e-4)):
        raise ValueError('Invalid predicted aligned error')


def write_json(path, data):
    partial = path.with_suffix('.partial')
    partial.write_text(json.dumps(data, indent=2) + '\n')
    partial.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--checkpoint', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--limit', type=int, default=128, help='Maximum new predictions in this production chunk')
    parser.add_argument('--max-length', type=int, default=512)
    parser.add_argument('--attention-chunk-size', type=int, default=16)
    args = parser.parse_args()
    if args.limit < 1 or args.max_length < 1 or args.attention_chunk_size < 1:
        raise ValueError('Positive chunk and length limits required')
    for directory in [args.inputs, args.checkpoint]:
        receipt = json.loads((directory / 'receipt.json').read_text())
        for name, checksum in receipt['artifacts'].items():
            if digest(directory / name) != checksum:
                raise ValueError('Changed source artifact: ' + name)
    if not torch.cuda.is_available():
        raise RuntimeError('An explicitly selected CUDA device is required')
    if not os.environ.get('CUDA_VISIBLE_DEVICES') or torch.cuda.device_count() != 1:
        raise RuntimeError('Select exactly one authorized GPU with CUDA_VISIBLE_DEVICES')
    torch.set_num_threads(4)
    torch.manual_seed(20260913)
    torch.cuda.manual_seed_all(20260913)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    args.output.mkdir(parents=True, exist_ok=True)
    lock = (args.output / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    config = {'script_sha256': digest(Path(__file__)),
        'input_receipt_sha256': digest(args.inputs / 'receipt.json'),
        'checkpoint_receipt_sha256': digest(args.checkpoint / 'receipt.json'),
        'python': platform.python_version(),
        'packages': {name: importlib.metadata.version(name) for name in ['torch', 'transformers', 'numpy', 'safetensors']},
        'modeling_source_sha256': digest(Path(modeling_esmfold.__file__)),
        'cuda_runtime': torch.version.cuda, 'gpu': torch.cuda.get_device_name(0),
        'visible_gpu': os.environ['CUDA_VISIBLE_DEVICES'], 'seed': 20260913,
        'dtype': 'ESM stem float16, folding trunk and heads float32; no autocast; TF32 disabled',
        'attention_chunk_size': args.attention_chunk_size, 'num_recycles_argument': None,
        'expected_trunk_passes': 4, 'batch_size': 1, 'max_length': args.max_length,
        'plddt_conversion': 'Transformers categorical_lddt returns 0..1; multiply by 100 for PDB and CA metadata',
        'mode': 'eval; torch.inference_mode; full single-chain sequence; no MSA; no template; no relaxation'}
    config_path = args.output / 'config.json'
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        raise ValueError('Existing prediction configuration differs; use a new output directory')
    write_json(config_path, config)
    sequences = list(read_fasta(args.inputs / 'candidates.faa'))
    if len({sid for sid, _ in sequences}) != len(sequences):
        raise ValueError('Duplicate prediction sequence IDs')
    eligible, deferred, cached = [], [], 0
    for sid, sequence in sequences:
        if sid != 'S' + hashlib.sha256(sequence.encode()).hexdigest():
            raise ValueError('Prediction sequence hash mismatch')
        if len(sequence) > args.max_length or not set(sequence) <= set(AA3.values()):
            deferred.append({'sequence_id': sid, 'length': len(sequence), 'status': 'length_deferred'
                             if len(sequence) > args.max_length else 'noncanonical_residue_deferred'})
            continue
        result_path = args.output / (sid + '.json')
        if result_path.exists():
            result = json.loads(result_path.read_text())
            if result['status'] != 'verified_prediction' or result['config_sha256'] != digest(config_path):
                raise ValueError('Invalid cached prediction receipt')
            for name, checksum in result['artifacts'].items():
                if digest(args.output / name) != checksum:
                    raise ValueError('Changed cached prediction artifact')
            cached += 1
        else:
            eligible.append((sid, sequence))
    write_json(args.output / 'deferred.json', deferred)
    print(json.dumps({'new_eligible': len(eligible), 'cached': cached, 'deferred': len(deferred),
                      'chunk_limit': args.limit}), flush=True)
    stop = [False]
    signal.signal(signal.SIGTERM, lambda *_: stop.__setitem__(0, True))
    start = time.monotonic()
    model = EsmForProteinFolding.from_pretrained(str(args.checkpoint), local_files_only=True,
                                               use_safetensors=True).eval()
    if model.config.esmfold_config.trunk.max_recycles != 4:
        raise ValueError('Unexpected checkpoint recycle configuration')
    model.esm.half()
    model.trunk.set_chunk_size(args.attention_chunk_size)
    model.cuda()
    print('model_loaded_seconds', round(time.monotonic() - start, 2), flush=True)
    completed, oom = 0, 0
    keys = ['positions', 'aatype', 'residue_index', 'plddt', 'atom37_atom_exists', 'residx_atom37_to_atom14']
    for sid, sequence in eligible[:args.limit]:
        if stop[0]:
            break
        torch.cuda.reset_peak_memory_stats()
        begin = time.monotonic()
        try:
            with torch.inference_mode():
                output = model.infer(sequence)
            torch.cuda.synchronize()
            seconds = time.monotonic() - begin
            pae = output['predicted_aligned_error'][0].float().cpu().numpy()
            raw_plddt = output['plddt'][0].float().cpu().numpy()
            if not np.isfinite(raw_plddt).all() or np.any(raw_plddt < 0) or np.any(raw_plddt > 1):
                raise ValueError('Unexpected native confidence scale')
            ca_plddt = raw_plddt[:, 1] * 100
            max_pae = float(output['max_predicted_aligned_error'].item())
            converted = {k: output[k] for k in keys}
            converted['plddt'] = converted['plddt'] * 100
            pdb = model.output_to_pdb(converted)[0]
            validate_prediction(sequence, pdb, pae, ca_plddt, max_pae)
            pdb_path = args.output / (sid + '.pdb')
            npz_path = args.output / (sid + '.npz')
            pdb_path.write_text(pdb)
            np.savez_compressed(npz_path, pae=pae, ca_plddt=ca_plddt,
                                max_pae=np.array(max_pae), sequence=np.array(sequence))
            result = {'status': 'verified_prediction', 'sequence_id': sid,
                'sequence_sha256': sid[1:], 'length': len(sequence), 'source': 'local ESMFold v1',
                'config_sha256': digest(config_path), 'inference_seconds': seconds,
                'peak_gpu_allocated_bytes': torch.cuda.max_memory_allocated(),
                'peak_gpu_reserved_bytes': torch.cuda.max_memory_reserved(),
                'mean_ca_plddt': float(ca_plddt.mean()),
                'fraction_ca_plddt_below50': float(np.mean(ca_plddt < 50)),
                'max_predicted_aligned_error': max_pae,
                'artifacts': {p.name: digest(p) for p in [pdb_path, npz_path]}}
            write_json(args.output / (sid + '.json'), result)
            completed += 1
            print(completed, sid, len(sequence), round(seconds, 2), round(result['mean_ca_plddt'], 2), flush=True)
            del output, converted
        except torch.cuda.OutOfMemoryError:
            oom += 1
            write_json(args.output / (sid + '.oom.json'), {'sequence_id': sid, 'length': len(sequence),
                       'status': 'oom_deferred', 'config_sha256': digest(config_path)})
            # Exit cleanly after recording the sequence; do not infer absence or truncate it.
            break
    write_json(args.output / 'last_chunk.json', {'status': 'production_chunk_finished',
        'new_predictions': completed, 'cached_predictions': cached, 'oom_deferred': oom,
        'remaining_eligible': len(eligible) - completed, 'length_or_alphabet_deferred': len(deferred),
        'interrupted': stop[0], 'config_sha256': digest(config_path)})


if __name__ == '__main__':
    main()
