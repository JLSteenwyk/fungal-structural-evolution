#!/usr/bin/env python3
"""Advance the completed reference core to OrthoFinder analysis of all 526 taxa."""
import csv
import fcntl
import json
import subprocess
import time
from pathlib import Path
from Bio import Phylo
from prepare_pfam import ROOT, digest


def taxa(path):
    names = [t.name for t in Phylo.read(path, 'newick').get_terminals()]
    if len(set(names)) != len(names):
        raise ValueError('Duplicate species-tree taxa')
    return set(names)


def main():
    control = ROOT / 'results/orthology/assignment-control-v2'
    control.mkdir(parents=True, exist_ok=True)
    lock = (control / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (control / 'run_config.json').exists():
        raise FileExistsError('Existing assignment requires explicit continuation review')
    core_receipt_path = ROOT / 'results/orthology/core-control-v1/receipt.json'
    core = json.loads(core_receipt_path.read_text())
    if core['returncode'] != 0 or core['status'] != 'core_finished_requires_full_assignment':
        raise ValueError('Successful core completion required')
    ir_path = ROOT / 'metadata/orthology_input_receipt.json'
    ir = json.loads(ir_path.read_text())
    manifest = ROOT / 'metadata/orthology_input_manifest.tsv'
    if digest(ir_path) != core['input_receipt_sha256'] or digest(manifest) != ir['input_manifest_sha256']:
        raise ValueError('Changed full-cohort inputs')
    rows = list(csv.DictReader(manifest.open(), delimiter='\t'))
    if len(rows) != 526 or len({r['taxon_id'] for r in rows}) != 526:
        raise ValueError('Full taxon scope changed')
    for row in rows:
        if digest(ROOT / row['input_path']) != row['source_sha256']:
            raise ValueError('Changed proteome: ' + row['taxon_id'])
    core_dir = ROOT / core['results_directory']
    core_tree = core_dir / 'Species_Tree/SpeciesTree_rooted.txt'
    if digest(core_tree) != core['species_tree_sha256'] or taxa(core_tree) != {r['taxon_id'] for r in rows if r['stage'] == 'core'}:
        raise ValueError('Core tree taxon/source mismatch')
    software_path = ROOT / 'metadata/orthofinder_software_receipt.json'
    if digest(software_path) != core['software_receipt_sha256']:
        raise ValueError('Changed OrthoFinder software receipt')
    software = json.loads(software_path.read_text())
    executable = ROOT / software['executable_path']
    if digest(executable) != software['executable_sha256']:
        raise ValueError('Changed OrthoFinder executable')
    dependency_lock = ROOT / 'environments/orthofinder-pip.lock.txt'
    runtime = subprocess.run([str(executable.parent / 'python'), '-c',
        'import numpy; assert numpy.__version__ == "2.2.6"; assert hasattr(numpy, "chararray"); print(numpy.__version__)'],
        capture_output=True, text=True, check=True).stdout.strip()
    protected = [core_tree, core_dir / 'Orthogroups/Orthogroups.tsv',
                 core_dir / 'Comparative_Genomics_Statistics/Statistics_Overall.tsv']
    core_hashes = {str(p.relative_to(ROOT)): digest(p) for p in protected}
    before = {p.name for p in core_dir.parent.iterdir() if p.is_dir()}
    command = [str(executable), '--assign', str(ROOT / 'data/orthology_inputs/v1/additional'),
               '--core', str(core_dir), '-n', 'full526v2', '-t', '32', '-a', '8',
               '-M', 'msa', '-S', 'diamond', '-A', 'famsa', '-T', 'fasttree']
    config = {'command': command, 'script_sha256': digest(Path(__file__)),
        'core_receipt_sha256': digest(core_receipt_path), 'core_artifacts': core_hashes,
        'input_receipt_sha256': digest(ir_path), 'software_receipt_sha256': digest(software_path),
        'dependency_lock_sha256': digest(dependency_lock), 'numpy_version': runtime,
        'expected_taxa': 526, 'additional_taxa': 462, 'results_parent': str(core_dir.parent.relative_to(ROOT)),
        'existing_result_directories': sorted(before),
        'scope': 'Assign all additional species, infer combined gene/species trees and reconciliation outputs; scientific validation remains separate.'}
    (control / 'run_config.json').write_text(json.dumps(config, indent=2) + '\n')
    begin = time.monotonic()
    with (control / 'stdout.log').open('w') as handle:
        completed = subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)
    after = [p for p in core_dir.parent.iterdir() if p.is_dir() and p.name not in before]
    trees = [p / 'Species_Tree/SpeciesTree_rooted.txt' for p in after
             if (p / 'Species_Tree/SpeciesTree_rooted.txt').exists()]
    result = {'returncode': completed.returncode, 'elapsed_seconds': time.monotonic() - begin,
              'config_sha256': digest(control / 'run_config.json'),
              'new_directories': [str(p.relative_to(ROOT)) for p in after],
              'status': 'failed_or_requires_output_review'}
    if any(digest(ROOT / p) != sha for p, sha in core_hashes.items()):
        raise ValueError('A protected core result changed during assignment')
    if completed.returncode == 0 and len(trees) == 1 and taxa(trees[0]) == {r['taxon_id'] for r in rows}:
        result.update(status='full_run_finished_requires_scientific_validation',
                      results_directory=str(trees[0].parents[1].relative_to(ROOT)),
                      species_tree_sha256=digest(trees[0]))
    (control / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)
    if result['status'] != 'full_run_finished_requires_scientific_validation':
        raise RuntimeError('OrthoFinder assignment failed; inspect log before any restart')


if __name__ == '__main__':
    main()
