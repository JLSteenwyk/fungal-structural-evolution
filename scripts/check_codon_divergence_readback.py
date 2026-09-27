#!/usr/bin/env python3
"""Ensure numerical and completeness checks reject rehashed comparison corruption."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    base = json.loads(Path('metadata/codon_divergence_identity_plan_20260927.json').read_text())
    for corruption in ('wrong_distance', 'missing_pair', 'wrong_omega'):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            copied = root / 'comparison'
            shutil.copytree(base['output'], copied)
            name = 'cases.tsv' if corruption == 'wrong_omega' else 'pairs.tsv'
            path = copied / name
            lines = path.read_text().splitlines()
            if corruption == 'missing_pair':
                del lines[1]
            else:
                fields = lines[1].split('\t')
                column = lines[0].split('\t').index('global_omega_delta' if corruption == 'wrong_omega' else 'original_ds_equal_alternative')
                fields[column] = str(float(fields[column]) + 1)
                lines[1] = '\t'.join(fields)
            path.write_text('\n'.join(lines)+'\n')
            plan = dict(base, output=str(copied))
            config = root / 'plan.json'
            config.write_text(json.dumps(plan))
            receipt_path = copied / 'receipt.json'
            receipt = json.loads(receipt_path.read_text())
            receipt['plan_sha256'] = sha(config)
            receipt['artifacts'][name] = sha(path)
            if corruption == 'missing_pair':
                receipt['pair_rows'] -= 1
            receipt_path.write_text(json.dumps(receipt))
            proof = root / 'proof.json'
            run = subprocess.run([sys.executable, 'scripts/readback_codon_alignment_divergence.py', '--plan', str(config), '--proof', str(proof)], capture_output=True, text=True)
            assert run.returncode != 0 and 'AssertionError' in run.stderr and 'check_rows' in run.stderr, run.stderr
            assert not proof.exists()
            print('Rejected', corruption, 'after updating artifact and plan hashes', flush=True)


if __name__ == '__main__':
    main()
