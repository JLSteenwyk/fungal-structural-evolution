#!/usr/bin/env python3
"""Run full expanded reconciliation guide sensitivities in isolated physical copies."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time
from audit_busco_gene_copies import ROOT, sha


def read(path):
    return json.loads(path.read_text())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args()
    plan, plan_sha = read(a.plan), sha(a.plan)
    output = ROOT / plan['output']
    if output.exists():
        raise FileExistsError('Inspect existing execution before recovery')

    def verify():
        if sha(a.plan) != plan_sha:
            raise ValueError('Execution plan changed')
        for name, digest in plan['pins'].items():
            if sha(ROOT / name) != digest:
                raise ValueError('Pinned source changed: ' + name)

    verify()
    inputs = ROOT / plan['inputs']
    receipt, audit = read(inputs / 'receipt.json'), read(inputs / 'readback.json')
    if (audit['status'] != 'passed_complete_expanded_input_readback'
            or audit['input_receipt_sha256'] != sha(inputs / 'receipt.json')
            or any(g['pending'] for g in audit['guides'])):
        raise ValueError('Complete independently validated inputs required')
    fixture = read(ROOT / plan['fixture_readback'])
    if fixture['status'] != 'passed_native_restart_fixture_output_readback':
        raise ValueError('Native command qualification required')
    output.mkdir(parents=True)
    state = dict(status='staging', plan_sha256=plan_sha, started_unix=time.time(), stages=[])

    def save():
        state['updated_unix'] = time.time()
        temporary = output / 'state.tmp'
        temporary.write_text(json.dumps(state, indent=2) + '\n')
        temporary.replace(output / 'state.json')

    save()
    try:
        for guide in receipt['guides']:
            verify()
            available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                                 if line.startswith('MemAvailable:'))) * 1024
            if available < plan['resources']['minimum_available_memory_gib'] * 2**30:
                raise RuntimeError('Insufficient currently available memory; inspect before retry')
            if shutil.disk_usage(ROOT).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
                raise RuntimeError('Insufficient disk headroom')
            name = guide['guide']
            state.update(status='staging', guide=name)
            save()
            old, base = inputs / name, output / name
            base.mkdir()
            manifest = old / 'copied_files.tsv'
            if sha(manifest) != guide['manifest_sha256']:
                raise ValueError('Changed input manifest')
            with manifest.open() as handle:
                copied = list(csv.DictReader(handle, delimiter='\t'))
            for row in copied:
                source, target = old / row['relative_path'], base / row['relative_path']
                if sha(source) != row['sha256']:
                    raise ValueError('Input snapshot changed')
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, target)
                if sha(target) != row['sha256'] or target.stat().st_ino == source.stat().st_ino:
                    raise ValueError('Fresh physical input copy failed')
            wd = base / 'Source/WorkingDirectory'
            descriptor = base / 'Source/Log.txt'
            descriptor.write_text('Full expanded conditional guide sensitivity; isolated execution inputs.\n'
                                  f'WorkingDirectory_Base: {wd}/\n'
                                  f'FN_Orthogroups: {wd}/clusters_OrthoFinder.txt_id_pairs.txt\n'
                                  f'WorkingDirectory_Trees: {wd}/\n')
            descriptor_sha = sha(descriptor)
            command = [str(ROOT / plan['executable'])] + plan.get('launcher_arguments', []) + ['--from-trees', str(base / 'Source'),
                       '-s', str(base / 'species_tree_taxa.nwk'), '-n', 'expanded_' + name,
                       '-t', str(plan['resources']['search_threads']),
                       '-a', str(plan['resources']['analysis_workers']), '-M', 'msa',
                       '-S', 'diamond', '-A', 'famsa', '-T', 'fasttree',
                       '--no-fix-files', '--save-space']
            env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', CUDA_VISIBLE_DEVICES='')
            started = time.time()
            with (base / 'stdout.log').open('w') as log:
                process = subprocess.Popen(command, cwd=base, env=env, stdout=log,
                                           stderr=subprocess.STDOUT, start_new_session=True)
                state.update(status='running_native_reconciliation', native_pid=process.pid,
                             command=command, native_started_unix=started)
                save()
                while process.poll() is None:
                    if shutil.disk_usage(ROOT).free < plan['resources']['emergency_free_disk_gib'] * 2**30:
                        os.killpg(process.pid, signal.SIGTERM)
                        try:
                            process.wait(timeout=30)
                        except subprocess.TimeoutExpired:
                            os.killpg(process.pid, signal.SIGKILL)
                            process.wait()
                        raise RuntimeError('Stopped native process group to preserve disk headroom')
                    time.sleep(20)
            if process.returncode != 0:
                raise RuntimeError('Native reconciliation exited ' + str(process.returncode))
            result = base / ('Results_expanded_' + name)
            required = ['Phylogenetic_Hierarchical_Orthogroups/N0.tsv',
                        'Gene_Duplication_Events/Duplications.tsv',
                        'Resolved_Gene_Trees/Resolved_Gene_Trees.txt',
                        'Species_Tree/SpeciesTree_rooted.txt',
                        'Comparative_Genomics_Statistics/Statistics_Overall.tsv']
            if any(not (result / f).is_file() or (result / f).stat().st_size == 0 for f in required):
                raise ValueError('Native exit-zero without required outputs')
            with (result / required[-1]).open() as handle:
                stats = {r[0]: r[1] for r in csv.reader(handle, delimiter='\t') if len(r) >= 2}
            if int(stats['Number of species']) != 526 or int(stats['Number of genes']) != 5815847:
                raise ValueError('Final statistics dimensions disagree')
            for row in copied:
                if sha(base / row['relative_path']) != row['sha256']:
                    raise ValueError('Native execution modified a copied input')
            if sha(descriptor) != descriptor_sha:
                raise ValueError('Native input descriptor changed')
            verify()
            stage = dict(guide=name, status='native_execution_complete_pending_full_output_readback',
                         elapsed_seconds=time.time() - started, result=str(result.relative_to(ROOT)),
                         mandatory_artifacts={f: sha(result / f) for f in required},
                         source_inputs_unchanged=True)
            state['stages'].append(stage)
            save()
        state['status'] = 'both_native_executions_complete_pending_full_output_readback'
        save()
    except Exception as exc:
        state.update(status='failed_requires_review', error=repr(exc))
        save()
        raise


if __name__ == '__main__':
    main()
