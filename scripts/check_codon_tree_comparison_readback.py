#!/usr/bin/env python3
"""Ensure independently reconstructed paths reject rehashed corrupt outputs."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    original = Path('results/cds/codon-tree-comparison-identity-check-20260927-v1')
    plan = json.loads(Path('metadata/codon_tree_comparison_identity_readback_plan_20260927.json').read_text())
    for corruption in ['wrong_distance', 'missing_pair']:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            copied = root / 'comparison'
            shutil.copytree(original, copied)
            path = copied / 'pairs.tsv'
            lines = path.read_text().splitlines()
            if corruption == 'wrong_distance':
                fields = lines[1].split('\t')
                fields[3] = str(float(fields[3]) + 1)
                lines[1] = '\t'.join(fields)
            else:
                del lines[1]
            path.write_text('\n'.join(lines) + '\n')
            receipt_path = copied / 'receipt.json'
            receipt = json.loads(receipt_path.read_text())
            receipt['artifacts']['pairs.tsv'] = hashlib.sha256(path.read_bytes()).hexdigest()
            if corruption == 'missing_pair':
                receipt['pair_rows'] -= 1
            receipt_path.write_text(json.dumps(receipt))
            plan['arguments']['output'] = str(copied)
            config = root / 'plan.json'
            config.write_text(json.dumps(plan))
            proof = root / 'proof.json'
            run = subprocess.run([sys.executable, 'scripts/readback_codon_tree_comparison.py',
                                  '--comparison-plan', str(config), '--output', str(proof)],
                                 capture_output=True, text=True)
            assert run.returncode != 0 and 'AssertionError' in run.stderr, run.stderr
            assert not proof.exists()
            print('Rejected', corruption, 'despite updated artifact hash/count', flush=True)


if __name__ == '__main__':
    main()
