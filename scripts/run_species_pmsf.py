#!/usr/bin/env python3
"""Run one supported full-matrix PMSF sensitivity using an audited frozen guide."""
import argparse
import fcntl
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
from Bio import Phylo, SeqIO
from audit_busco_gene_copies import ROOT, sha


def checked(folder):
    row = json.loads((folder / 'receipt.json').read_text())
    for name, digest in row['artifacts'].items():
        if sha(folder / name) != digest:
            raise ValueError('Changed input artifact: ' + name)
    return row


def validate_tree(path, expected):
    tree = Phylo.read(path, 'newick')
    tips = [x.name for x in tree.get_terminals()]
    if len(tips) != len(set(tips)) or set(tips) != expected:
        raise ValueError('PMSF tree tips differ')
    if any(n.branch_length is not None and (not math.isfinite(n.branch_length) or n.branch_length < 0)
           for n in tree.find_clades()):
        raise ValueError('Invalid PMSF branch length')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['matrix','guide-audit','resources','output']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Inspect existing execution before recovery; use a new output')
    lock = (ROOT / 'results/phylogeny/.species_pmsf.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    matrix = checked(a.matrix); guide = checked(a.guide_audit); resources = checked(a.resources)
    if guide['status'] != 'passed_full_species_guide_readback' or guide['taxa'] != 526:
        raise ValueError('Full audited guide required')
    if resources['status'] != 'prepared_full_matrix_mixture_resource_scenarios':
        raise ValueError('Resource assessment required')
    if not any(row['matrix_receipt_sha256'] == sha(a.matrix / 'receipt.json') for row in resources['matrices']):
        raise ValueError('Matrix outside resource assessment')
    with (a.matrix / 'matrix.faa').open() as handle:
        records = list(SeqIO.parse(handle, 'fasta'))
    taxa = {r.id for r in records}
    if len(records) != 526 or len(taxa) != 526 or any(len(r.seq) != matrix['columns'] for r in records):
        raise ValueError('Full matrix dimensions differ')
    validate_tree(a.guide_audit / 'guide.treefile', taxa)
    plan = resources['proposed_design']
    available_kib = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))
    if available_kib < plan['minimum_available_memory_gib'] * 1024**2:
        raise RuntimeError('Insufficient currently available memory; no launch')
    if shutil.disk_usage(ROOT).free < plan['output_allowance_gb'] * 10**9:
        raise RuntimeError('Insufficient output headroom; no launch')
    executable = Path(json.loads((a.guide_audit / 'run_config.json').read_text())['command'][0])
    if sha(executable) != '40424ccdb1d79c304641f910cb6c172ebb670214e50958352018e3ff9906ab8f':
        raise ValueError('IQ-TREE build differs')
    a.output.mkdir(parents=True)
    frozen_guide = a.output / 'input_guide.treefile'
    shutil.copyfile(a.guide_audit / 'guide.treefile', frozen_guide)
    command = [str(executable), '-s', str((a.matrix / 'matrix.faa').resolve()), '-st', 'AA',
        '-m', plan['model'], '--tree-freq', str(frozen_guide.resolve()),
        '-T', str(plan['threads']), '--mem', plan['iqtree_memory_limit'], '--seed', '20260913',
        '--alrt', '1000', '-B', '1000', '--bnni', '--boot-trees',
        '--prefix', str((a.output / 'pmsf').resolve())]
    pins = {str(path.resolve()): sha(path) for path in [Path(__file__),
        Path(__file__).with_name('audit_busco_gene_copies.py'), a.matrix / 'receipt.json',
        a.matrix / 'matrix.faa', a.guide_audit / 'receipt.json', a.guide_audit / 'guide.treefile',
        a.resources / 'receipt.json', executable, frozen_guide]}
    config = {'command': command, 'pinned_files': pins, 'resource_plan': plan,
        'available_memory_kib_at_launch': available_kib,
        'matrix_receipt_sha256': sha(a.matrix / 'receipt.json'),
        'guide_audit_receipt_sha256': sha(a.guide_audit / 'receipt.json'),
        'interpretation': 'One full-taxon alignment/guide combination in crossed PMSF sensitivity; support conditions on estimated profiles. No final species-tree, adequacy or dating claim.'}
    config_path = a.output / 'config.json'; config_path.write_text(json.dumps(config,indent=2)+'\n')
    print('Launching supported full-matrix PMSF', matrix['columns'], 'sites', flush=True)
    start = time.monotonic(); env = os.environ.copy(); env['OPENBLAS_NUM_THREADS']='1'
    with (a.output / 'stdout.log').open('w') as log:
        run = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    result = {'status': 'pmsf_execution_requires_review', 'returncode': run.returncode,
              'elapsed_seconds': time.monotonic()-start, 'config_sha256': sha(config_path)}
    if run.returncode == 0:
        for path, expected in pins.items():
            if sha(Path(path)) != expected:
                raise ValueError('Pinned input changed during inference')
        validate_tree(a.output / 'pmsf.treefile', taxa)
        n = 0
        with (a.output / 'pmsf.ufboot').open() as handle:
            for tree in Phylo.parse(handle, 'newick'):
                names=[tip.name for tip in tree.get_terminals()]
                if len(names) != 526 or set(names) != taxa:
                    raise ValueError('Bootstrap taxon grid differs')
                n += 1
        if n != 1000:
            raise ValueError('Bootstrap replicate count differs')
        if not (a.output / 'pmsf.sitefreq').is_file():
            raise ValueError('PMSF profiles missing')
        result.update(status='complete_pmsf_execution_pending_full_audit', taxa=526,
            columns=matrix['columns'], bootstrap_trees=n,
            artifacts={p.name:sha(p) for p in a.output.iterdir() if p.is_file()},
            interpretation='Execution and tip grids checked; full support/profile/model-output audit remains required. This is one sensitivity analysis, not the final species phylogeny.')
    (a.output / 'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    if run.returncode:
        raise SystemExit(run.returncode)


if __name__ == '__main__':
    main()
