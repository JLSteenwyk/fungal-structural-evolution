#!/usr/bin/env python3
"""Verify checkpoint loading and ensure the unused contact head is absent from folding execution."""
import hashlib
import json
import os
from pathlib import Path
import numpy as np
import torch
from transformers import EsmForProteinFolding
from transformers.models.esm import modeling_esmfold
from prepare_pfam import ROOT, digest


def main():
    checkpoint = ROOT / 'data/prediction_models/esmfold-v1'
    predictions = ROOT / 'results/predictions/esmfold-marker-v1'
    receipt = json.loads((checkpoint / 'receipt.json').read_text())
    for name, expected in receipt['artifacts'].items():
        if digest(checkpoint / name) != expected:
            raise ValueError('Changed checkpoint')
    if not os.environ.get('CUDA_VISIBLE_DEVICES') or torch.cuda.device_count() != 1:
        raise RuntimeError('Explicitly select one available GPU')
    torch.set_num_threads(4)
    torch.manual_seed(20260913)
    torch.cuda.manual_seed_all(20260913)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    model, loading = EsmForProteinFolding.from_pretrained(str(checkpoint), local_files_only=True,
                         use_safetensors=True, output_loading_info=True)
    expected = {'esm.contact_head.regression.bias', 'esm.contact_head.regression.weight'}
    if (set(loading['missing_keys']) != expected or loading['unexpected_keys']
            or loading.get('mismatched_keys') or loading.get('error_msgs')):
        raise ValueError('Unexpected checkpoint loading discrepancy: ' + repr(loading))
    def reject_contact_head(*_):
        raise RuntimeError('Uninitialized contact head was invoked during folding')
    hook = model.esm.contact_head.register_forward_pre_hook(reject_contact_head)
    model.eval()
    model.esm.half()
    model.trunk.set_chunk_size(64)
    model.cuda()
    path = next(p for p in sorted(predictions.glob('S*.json')) if not p.name.endswith('.oom.json'))
    row = json.loads(path.read_text())
    for name, expected_hash in row['artifacts'].items():
        if digest(predictions / name) != expected_hash:
            raise ValueError('Changed production artifact')
    with np.load(predictions / (row['sequence_id'] + '.npz'), allow_pickle=False) as stored:
        sequence = str(stored['sequence'])
        with torch.inference_mode():
            output = model.infer(sequence)
        pae = output['predicted_aligned_error'][0].float().cpu().numpy()
        confidence = output['plddt'][0, :, 1].float().cpu().numpy() * 100
        pae_difference = float(np.max(np.abs(pae - stored['pae'])))
        confidence_difference = float(np.max(np.abs(confidence - stored['ca_plddt'])))
        if not np.allclose(pae, stored['pae'], rtol=0, atol=1e-4) or not np.allclose(confidence, stored['ca_plddt'], rtol=0, atol=.01):
            raise ValueError('Repeated inference differs beyond declared tolerance')
    hook.remove()
    source_dir = Path(modeling_esmfold.__file__).parent
    result = {'status': 'loading_and_contact_head_execution_audit_passed',
        'loading_info': loading, 'contact_head_forward_hook': 'raises if invoked; complete folding inference passed',
        'repeated_sequence_id': row['sequence_id'], 'length': len(sequence),
        'max_absolute_pae_difference': pae_difference, 'max_absolute_ca_plddt_difference': confidence_difference,
        'production_prediction_receipt_sha256': digest(path),
        'production_config_sha256': digest(predictions / 'config.json'),
        'checkpoint_receipt_sha256': digest(checkpoint / 'receipt.json'),
        'script_sha256': digest(Path(__file__)),
        'esm_source_files': {str(p.relative_to(source_dir)): digest(p) for p in sorted(source_dir.rglob('*.py'))},
        'interpretation': 'Unused contact-regression parameters do not participate in this executed folding path. Repeatability of one sequence does not establish experimental accuracy or portability.'}
    (ROOT / 'metadata/esmfold_execution_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], pae_difference, confidence_difference)


if __name__ == '__main__':
    main()
