#!/usr/bin/env python3
"""Run a fresh, isolated MAFFT-guide expanded-reconciliation recovery."""
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = read(args.plan)
    plan_sha = sha(args.plan)
    output = ROOT / plan['output']
    if output.exists():
        raise FileExistsError('Recovery output already exists: ' + str(output))

    def verify_pins():
        if sha(args.plan) != plan_sha:
            raise ValueError('Execution plan changed during recovery')
        for name, digest in plan['pins'].items():
            if sha(ROOT / name) != digest:
                raise ValueError('Pinned source changed: ' + name)

    def available_bytes():
        for line in Path('/proc/meminfo').read_text().splitlines():
            if line.startswith('MemAvailable:'):
                return int(line.split()[1]) * 1024
        raise RuntimeError('MemAvailable unavailable')

    verify_pins()
    if available_bytes() < plan['resources']['minimum_available_memory_gib'] * 2**30:
        raise RuntimeError('Insufficient available memory before recovery')
    if shutil.disk_usage(ROOT).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
        raise RuntimeError('Insufficient disk headroom before recovery')

    inputs = ROOT / plan['inputs']
    input_receipt = read(inputs / 'receipt.json')
    input_readback = read(inputs / 'readback.json')
    if (input_readback['status'] != 'passed_complete_expanded_input_readback'
            or input_readback['input_receipt_sha256'] != sha(inputs / 'receipt.json')):
        raise ValueError('Complete independent input readback is required')
    guide = next((row for row in input_receipt['guides'] if row['guide'] == 'mafft'), None)
    checked = next((row for row in input_readback['guides'] if row['guide'] == 'mafft'), None)
    if guide is None or checked is None or checked['pending']:
        raise ValueError('Complete MAFFT-guide inputs are required')
    fixture = read(ROOT / plan['fixture_readback'])
    if fixture['status'] != 'passed_native_restart_fixture_output_readback':
        raise ValueError('Native restart qualification is required')

    source = inputs / 'mafft'
    manifest = source / 'copied_files.tsv'
    if sha(manifest) != guide['manifest_sha256']:
        raise ValueError('MAFFT input manifest changed')
    with manifest.open() as handle:
        copied = list(csv.DictReader(handle, delimiter='\t'))

    output.mkdir(parents=True)
    state = {
        'status': 'staging', 'guide': 'mafft', 'plan_sha256': plan_sha,
        'started_unix': time.time(), 'input_manifest_sha256': sha(manifest),
        'input_files': len(copied),
    }

    def save_state():
        state['updated_unix'] = time.time()
        temporary = output / 'state.tmp'
        temporary.write_text(json.dumps(state, indent=2) + '\n')
        temporary.replace(output / 'state.json')

    save_state()
    try:
        for row in copied:
            original = source / row['relative_path']
            target = output / row['relative_path']
            if sha(original) != row['sha256']:
                raise ValueError('Changed input file: ' + row['relative_path'])
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(original, target)
            if sha(target) != row['sha256'] or target.stat().st_ino == original.stat().st_ino:
                raise ValueError('Fresh physical input copy failed: ' + row['relative_path'])

        workdir = output / 'Source/WorkingDirectory'
        descriptor = output / 'Source/Log.txt'
        descriptor.write_text(
            'Fresh isolated MAFFT-guide expanded reconciliation recovery.\n'
            f'WorkingDirectory_Base: {workdir}/\n'
            f'FN_Orthogroups: {workdir}/clusters_OrthoFinder.txt_id_pairs.txt\n'
            f'WorkingDirectory_Trees: {workdir}/\n'
        )
        descriptor_sha = sha(descriptor)
        command = [
            str(ROOT / plan['executable']), str(ROOT / 'scripts/run_orthofinder_with_progress.py'),
            '--stall-timeout-seconds', str(plan['task_stall_allowance_seconds']), '--',
            '--from-trees', str(output / 'Source'), '-s', str(output / 'species_tree_taxa.nwk'),
            '-n', plan['result_name'], '-t', str(plan['resources']['search_threads']),
            '-a', str(plan['resources']['analysis_workers']), '-M', 'msa', '-S', 'diamond',
            '-A', 'famsa', '-T', 'fasttree', '--no-fix-files', '--save-space',
        ]
        environment = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
                           MKL_NUM_THREADS='1', CUDA_VISIBLE_DEVICES='')
        started = time.time()
        with (output / 'stdout.log').open('w') as log:
            process = subprocess.Popen(command, cwd=output, env=environment, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            state.update(status='running_native_reconciliation', native_pid=process.pid,
                         native_started_unix=started, command=command)
            save_state()
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

        result = output / ('Results_' + plan['result_name'])
        required = [
            'Phylogenetic_Hierarchical_Orthogroups/N0.tsv',
            'Gene_Duplication_Events/Duplications.tsv',
            'Resolved_Gene_Trees/Resolved_Gene_Trees.txt',
            'Species_Tree/SpeciesTree_rooted.txt',
            'Comparative_Genomics_Statistics/Statistics_Overall.tsv',
        ]
        if any(not (result / name).is_file() or (result / name).stat().st_size == 0 for name in required):
            raise ValueError('Native exit-zero without required output')
        with (result / required[-1]).open() as handle:
            statistics = {row[0]: row[1] for row in csv.reader(handle, delimiter='\t') if len(row) >= 2}
        if int(statistics['Number of species']) != 526 or int(statistics['Number of genes']) != 5815847:
            raise ValueError('Final native dimensions disagree with frozen panel')
        for row in copied:
            if sha(output / row['relative_path']) != row['sha256']:
                raise ValueError('Native execution changed copied input: ' + row['relative_path'])
        if sha(descriptor) != descriptor_sha:
            raise ValueError('Native descriptor changed')
        verify_pins()
        receipt = {
            'status': 'complete_native_mafft_reconciliation_pending_independent_output_readback',
            'guide': 'mafft', 'plan_sha256': plan_sha, 'input_manifest_sha256': sha(manifest),
            'input_files': len(copied), 'elapsed_seconds': time.time() - started,
            'result': str(result.relative_to(ROOT)), 'source_inputs_unchanged': True,
            'mandatory_artifacts': {name: sha(result / name) for name in required},
            'statistics': {'species': int(statistics['Number of species']), 'genes': int(statistics['Number of genes'])},
            'interpretation': 'Fresh isolated native MAFFT-guide execution. Independent output readback and guide-sensitivity reconciliation review remain required.',
        }
        (output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        state.update(status=receipt['status'], native_returncode=process.returncode,
                     receipt_sha256=sha(output / 'receipt.json'))
        save_state()
    except Exception as exc:
        state.update(status='failed_requires_review', error=repr(exc))
        save_state()
        raise


if __name__ == '__main__':
    main()
