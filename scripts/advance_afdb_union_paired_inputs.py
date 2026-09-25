#!/usr/bin/env python3
"""Prepare and independently read back paired inputs from an audited AFDB union."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from compare_marker_structures import sha
from assess_pae_sensitivity import checked_receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args(); plan = json.loads(a.plan.read_text()); digest = sha(a.plan)
    def verify():
        if sha(a.plan) != digest:
            raise ValueError('Plan changed')
        for path, expected in plan['pins'].items():
            if sha(Path(path)) != expected:
                raise ValueError('Changed dependency: ' + path)
    verify()
    audit = json.loads(Path(plan['union_readback']).read_text())
    if audit['status'] != 'passed_full_disjoint_afdb_union_row_readback' or audit['models'] != plan['models']:
        raise ValueError('Union readback missing or incomplete')
    for path, expected in audit['source_sha256'].items():
        if sha(Path(path)) != expected:
            raise ValueError('Changed union audit source')
    mapping, encodings = Path(plan['mapping']), Path(plan['encodings'])
    mr, er = checked_receipt(mapping), checked_receipt(encodings)
    if er['mapping_receipt_sha256'] != sha(mapping / 'receipt.json') or er['models'] != plan['models']:
        raise ValueError('Union mapping and encoding scope differ')
    if er['confidence_stage'] != 'plddt_and_pae':
        raise ValueError('Unqualified union')
    control = Path(plan['control']); inputs = Path(plan['inputs']); readback = Path(plan['readback'])
    if any(x.exists() for x in [control, inputs, readback]):
        raise FileExistsError('Use fresh outputs; inspect existing run before recovery')
    control.mkdir(parents=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', CUDA_VISIBLE_DEVICES='')
    commands = [('preparation', ['scripts/prepare_paired_phylogenetic_inputs.py', '--encodings', str(encodings), '--snapshot', str(mapping), '--matrix', plan['matrix'], '--output', str(inputs)]),
                ('readback', ['scripts/readback_paired_inputs_from_encodings.py', '--inputs', str(inputs), '--output', str(readback)])]
    for stage, command in commands:
        verify()
        if shutil.disk_usage(control).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
            raise ValueError('Insufficient disk')
        (control / 'state.json').write_text(json.dumps(dict(status='running', stage=stage, plan_sha256=digest)) + '\n')
        with (control / (stage + '.log')).open('w') as handle:
            subprocess.run([sys.executable] + command, env=env, stdout=handle, stderr=subprocess.STDOUT, check=True)
    result = checked_receipt(inputs); read = json.loads(readback.read_text()); baseline = checked_receipt(Path(plan['baseline']))
    if read['status'] != 'passed_complete_paired_inputs_from_qualified_arrays_readback' or read['source_receipt_sha256'] != sha(inputs / 'receipt.json'):
        raise ValueError('Readback binding differs')
    for field in ['mask', 'eligibility']:
        if result[field] != baseline[field]:
            raise ValueError('Changed paired filtering')
    if result['source_receipts']['matrix']['sha256'] != baseline['source_receipts']['matrix']['sha256']:
        raise ValueError('Changed matrix')
    verify()
    receipt = dict(status='complete_afdb_union_paired_inputs_and_readback', plan_sha256=digest,
                   paired_receipt_sha256=sha(inputs / 'receipt.json'), readback_sha256=sha(readback),
                   ready_markers=result['ready_markers'], scope='Full union paired inputs and character readback; existing evolutionary fits are not updated by this preparation.')
    (control / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    (control / 'state.json').write_text(json.dumps(receipt) + '\n')


if __name__ == '__main__':
    main()
