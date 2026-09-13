#!/usr/bin/env python3
"""Run the reference-core dependency of full-cohort OrthoFinder inference."""
import csv
import fcntl
import hashlib
import json
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    control = ROOT / 'results/orthology/core-control-v1'
    control.mkdir(parents=True, exist_ok=True)
    lock = (control / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    inputs = json.loads((ROOT / 'metadata/orthology_input_receipt.json').read_text())
    manifest = ROOT / 'metadata/orthology_input_manifest.tsv'
    if sha(manifest) != inputs['input_manifest_sha256']:
        raise ValueError('Changed orthology input manifest')
    with manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len(rows) != inputs['taxa'] or sum(r['stage'] == 'core' for r in rows) != inputs['core_taxa']:
        raise ValueError('Core/full-cohort scope mismatch')
    for row in rows:
        if sha(ROOT / row['input_path']) != row['source_sha256']:
            raise ValueError(f'Changed input: {row["taxon_id"]}')
    software = json.loads((ROOT / 'metadata/orthofinder_software_receipt.json').read_text())
    executable = ROOT / software['executable_path']
    if sha(executable) != software['executable_sha256']:
        raise ValueError('Changed OrthoFinder entrypoint')
    destination = ROOT / 'results/orthology/core-v1'
    if destination.exists():
        raise FileExistsError('Existing core output requires explicit restart review; never force overwrite')
    command = [str(executable), '-f', str(ROOT / 'data/orthology_inputs/v1/core'),
               '-o', str(destination), '-t', '32', '-a', '8', '-M', 'msa',
               '-S', 'diamond', '-A', 'famsa', '-T', 'fasttree']
    config = {'command': command, 'input_receipt_sha256': sha(ROOT / 'metadata/orthology_input_receipt.json'),
              'software_receipt_sha256': sha(ROOT / 'metadata/orthofinder_software_receipt.json'),
              'scope': '64 reference taxa; all remaining 462 taxa require subsequent assignment and combined analysis'}
    (control / 'run_config.json').write_text(json.dumps(config, indent=2) + '\n')
    start = time.monotonic()
    with (control / 'stdout.log').open('w') as log:
        result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    receipt = dict(config, returncode=result.returncode, elapsed_seconds=time.monotonic() - start)
    trees = list(destination.glob('**/Species_Tree/SpeciesTree_rooted.txt'))
    if result.returncode == 0 and len(trees) == 1:
        receipt.update(status='core_finished_requires_full_assignment', results_directory=str(trees[0].parents[1].relative_to(ROOT)),
                       species_tree_sha256=sha(trees[0]))
    else:
        receipt['status'] = 'failed_or_requires_output_review'
    (control / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(receipt['status'], flush=True)
    if result.returncode:
        raise RuntimeError(f'OrthoFinder exited {result.returncode}')


if __name__ == '__main__':
    main()
