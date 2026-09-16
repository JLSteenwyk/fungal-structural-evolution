#!/usr/bin/env python3
"""Exercise installed OrthoFinder restart on isolated, explicitly synthetic data."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def snapshot(folder):
    return {str(p.relative_to(folder)): sha(p) for p in folder.rglob('*') if p.is_file()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    executable = ROOT / '.cache/envs/orthofinder/bin/orthofinder'
    cases = []
    for case in ['missing_ids_tree', 'prepared_ids_tree']:
        folder = out / case
        source = folder / 'Source'
        wd = source / 'WorkingDirectory'
        trees = wd / 'Trees_ids'
        trees.mkdir(parents=True)
        (wd / 'SpeciesIDs.txt').write_text(''.join(f'{i}: Taxon{i}.faa\n' for i in range(4)))
        (wd / 'SequenceIDs.txt').write_text(''.join(f'{i}_{g}: protein_{i}_{g}\n' for i in range(4) for g in range(4 if i == 0 else 3)))
        for i in range(4):
            (wd / f'Species{i}.fa').write_text(''.join(f'>{i}_{g}\n' + 'ACDEFGHIKLMNPQRSTVWY' * 3 + '\n' for g in range(4 if i == 0 else 3)))
        clusters = wd / 'clusters_OrthoFinder.txt_id_pairs.txt'
        clusters.write_text('(mclmatrixbegin\n' + ''.join(str(g) + ' ' + ' '.join(f'{i}_{g}' for i in range(4)) + (' 0_3' if g == 2 else '') + ' $\n' for g in range(3)) + ')\n')
        for g in range(3):
            left = f'0_{g}:0.1' if g != 2 else '(0_2:0.05,0_3:0.05)1:0.05'
            (trees / f'OG{g:07d}.txt').write_text(f'(({left},1_{g}:0.1)1:0.1,(2_{g}:0.1,3_{g}:0.1)1:0.1);\n')
        (source / 'Log.txt').write_text('Synthetic native restart fixture; not a scientific dataset.\n' +
            f'WorkingDirectory_Base: {wd}/\nFN_Orthogroups: {clusters}\nWorkingDirectory_Trees: {wd}/\n')
        user_tree = folder / 'species_tree.nwk'
        user_tree.write_text('((Taxon0:0.1,Taxon1:0.1):0.1,(Taxon2:0.1,Taxon3:0.1):0.1);\n')
        if case == 'prepared_ids_tree':
            (wd / 'SpeciesTree_unrooted_ids.txt').write_text('((0:0.1,1:0.1):0.1,(2:0.1,3:0.1):0.1);\n')
        before = snapshot(source)
        command = [str(executable), '--from-trees', str(source), '-s', str(user_tree),
                   '-n', 'native_fixture', '-t', '2', '-a', '1', '-M', 'msa',
                   '-S', 'diamond', '-A', 'famsa', '-T', 'fasttree', '--no-fix-files']
        start = time.time()
        timed_out = False
        with (folder / 'stdout.log').open('w') as log:
            p = subprocess.Popen(command, cwd=folder, stdout=log, stderr=subprocess.STDOUT,
                                 start_new_session=True,
                                 env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1'))
            try:
                rc = p.wait(timeout=120)
            except subprocess.TimeoutExpired:
                timed_out = True
                os.killpg(p.pid, signal.SIGTERM)
                try:
                    rc = p.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid, signal.SIGKILL)
                    rc = p.wait()
        after = snapshot(source)
        result = dict(case=case, returncode=rc, timed_out=timed_out,
                      elapsed_seconds=time.time() - start, command=command,
                      source_files_changed=[name for name, h in before.items() if after.get(name) != h],
                      source_files_added=sorted(set(after) - set(before)),
                      source_before=before, source_after=after,
                      log_sha256=sha(folder / 'stdout.log'),
                      new_result_directories=[str(p.relative_to(out)) for p in folder.glob('Results_*')])
        cases.append(result)
        print(case, rc, 'changed', result['source_files_changed'], 'added', result['source_files_added'], flush=True)
    receipt = dict(status='completed_native_restart_behavior_observation', cases=cases,
                   script_sha256=sha(Path(__file__)), executable_sha256=sha(executable),
                   interpretation='Two disposable four-taxon/three-family/13-protein fixtures using the actual installed CLI. Tests restart plumbing and observed source writes, not scientific model validity or full-dataset reconciliation.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
