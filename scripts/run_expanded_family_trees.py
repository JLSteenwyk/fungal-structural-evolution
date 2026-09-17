#!/usr/bin/env python3
"""Run isolated, checkpointed FAMSA/native-trim/FastTree tasks for new families."""
import argparse
import concurrent.futures as futures
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import time
from Bio import SeqIO, Phylo
from orthofinder.tools import trim


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def save(path, value):
    tmp = path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def sequences(path):
    records = list(SeqIO.parse(path, 'fasta'))
    result = {r.id: str(r.seq).upper() for r in records}
    if len(result) != len(records) or not records:
        raise ValueError('Empty or duplicate FASTA identifiers')
    return result


def check_alignment(source, raw, aligned):
    original, before, after = sequences(source), sequences(raw), sequences(aligned)
    if set(original) != set(before) or set(before) != set(after):
        raise ValueError('Alignment identifier loss or addition')
    if len({len(s) for s in before.values()}) != 1 or len({len(s) for s in after.values()}) != 1:
        raise ValueError('Nonrectangular alignment')
    if any(before[k].replace('-', '') != original[k] for k in original):
        raise ValueError('Alignment changed source residues')
    names = sorted(original)
    # Native trimming treats stop symbols as gaps when it rewrites columns.
    normalize = lambda s: s.replace('*', '-')
    raw_columns = iter(zip(*(normalize(before[k]) for k in names)))
    for column in zip(*(normalize(after[k]) for k in names)):
        for candidate in raw_columns:
            if candidate == column:
                break
        else:
            raise ValueError('Trimmed columns not an ordered subset of raw columns')
    if any(not s.replace('-', '').replace('*', '') for s in after.values()):
        raise ValueError('Trimming left an empty sequence')
    return dict(proteins=len(original), raw_columns=len(next(iter(before.values()))),
                trimmed_columns=len(next(iter(after.values()))),
                raw_residues=sum(len(s.replace('-', '').replace('*', '')) for s in before.values()),
                trimmed_residues=sum(len(s.replace('-', '').replace('*', '')) for s in after.values()))


def check_tree(path, source):
    text = path.read_text().strip()
    if text.count(';') != 1 or not text.endswith(';'):
        raise ValueError('Incomplete or multiple trees')
    tree = Phylo.read(path, 'newick')
    tips = [n.name for n in tree.get_terminals()]
    if len(tips) != len(set(tips)) or set(tips) != set(sequences(source)):
        raise ValueError('Tree tip universe differs')
    for node in tree.find_clades():
        if node is tree.root:
            continue
        if node.branch_length is None or not math.isfinite(node.branch_length) or node.branch_length < 0:
            raise ValueError('Missing, negative or nonfinite branch length')
    return len(tips)


def command(argv, log_path, timeout, stdout_path=None):
    started = time.time()
    with log_path.open('w') as log:
        output = stdout_path.open('w') if stdout_path else log
        try:
            process = subprocess.Popen(argv, stdout=output, stderr=log, start_new_session=True)
            try:
                rc = process.wait(timeout=timeout)
            except BaseException:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                raise
            if rc != 0:
                raise RuntimeError(f'Command exited {rc}: {argv}')
        finally:
            if stdout_path:
                output.close()
    return dict(argv=argv, seconds=time.time() - started)


def run_task(row, plan, plan_hash):
    key = row['membership_sha256']
    folder = Path(plan['output']) / 'families' / key
    source = Path(plan['inputs']) / row['relative_path']
    if sha(source) != row['sha256']:
        raise ValueError('Input FASTA changed')
    receipt = folder / 'receipt.json'
    if receipt.exists():
        old = json.loads(receipt.read_text())
        if old['plan_sha256'] != plan_hash or old['source_sha256'] != row['sha256']:
            raise ValueError('Checkpoint provenance differs')
        for name, digest in old['artifacts'].items():
            if sha(folder / name) != digest:
                raise ValueError('Checkpoint artifact changed')
        return old
    if folder.exists():
        raise ValueError('Partial task requires review: ' + key)
    folder.mkdir(parents=True)
    raw, aligned, tree = (folder / name for name in ['raw.faa', 'trimmed.faa', 'tree.nwk'])
    timeout = plan['resources']['step_timeout_hours'] * 3600
    calls = [command([plan['famsa'], str(source), str(raw), '-t', '1'], folder / 'famsa.log', timeout)]
    # Separate output retains the untrimmed alignment for independent checks.
    trim.main(str(raw), str(aligned), 10, 0.1, 500, 0.75, False)
    dimensions = check_alignment(source, raw, aligned)
    calls.append(command([plan['fasttree'], str(aligned)], folder / 'fasttree.log', timeout, tree))
    tips = check_tree(tree, source)
    if tips != int(row['proteins']) or sha(source) != row['sha256']:
        raise ValueError('Input changed during inference')
    result = dict(status='complete_family_alignment_tree_requires_native_readback',
                  membership_sha256=key, plan_sha256=plan_hash, source_sha256=row['sha256'],
                  dimensions=dimensions, commands=calls,
                  artifacts={p.name: sha(p) for p in folder.iterdir() if p.is_file()})
    save(receipt, result)
    return result


def worker_init(memory_gib):
    resource.setrlimit(resource.RLIMIT_AS, (memory_gib * 2**30, memory_gib * 2**30))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    a = ap.parse_args()
    plan = json.loads(a.plan.read_text())
    for path, digest in plan['pins'].items():
        if sha(path) != digest:
            raise ValueError('Changed pinned file: ' + path)
    root = Path(plan['output'])
    root.mkdir(parents=True, exist_ok=True)
    lock = (root / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    plan_hash = sha(a.plan)
    stamp = root / 'plan_sha256.txt'
    if stamp.exists() and stamp.read_text().strip() != plan_hash:
        raise ValueError('Existing output has another plan')
    stamp.write_text(plan_hash + '\n')
    inputs = Path(plan['inputs'])
    receipt = json.loads((inputs / 'receipt.json').read_text())
    audit = json.loads((inputs / 'readback.json').read_text())
    if audit['status'] != 'passed_complete_new_family_sequence_readback' or audit['staging_receipt_sha256'] != sha(inputs / 'receipt.json'):
        raise ValueError('Independent sequence readback required')
    if sha(inputs / 'family_sequences.tsv') != receipt['artifacts']['family_sequences.tsv']:
        raise ValueError('Changed sequence table')
    with (inputs / 'family_sequences.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len(rows) != plan['expected_families'] or len({r['membership_sha256'] for r in rows}) != len(rows):
        raise ValueError('Task universe differs')
    rows.sort(key=lambda r: (-int(r['residues']), r['membership_sha256']))
    limits = plan['resources']
    state = dict(status='running', pid=os.getpid(), plan_sha256=plan_hash,
                 total=len(rows), completed=0, failures=[], started_at=time.time())
    worker_count = limits['workers']
    iterator = iter(rows)
    pending = {}
    with futures.ProcessPoolExecutor(max_workers=worker_count, initializer=worker_init,
                                     initargs=(limits['per_worker_address_space_gib'],)) as pool:
        exhausted = False
        while pending or not exhausted:
            while len(pending) < worker_count and not exhausted:
                if shutil.disk_usage(root).free < limits['minimum_free_disk_gib'] * 2**30:
                    raise RuntimeError('Free disk gate')
                with Path('/proc/meminfo').open() as handle:
                    available = next(int(line.split()[1]) * 1024 for line in handle if line.startswith('MemAvailable:'))
                if available < limits['minimum_available_memory_gib'] * 2**30:
                    # Live jobs finish; new jobs wait for available memory.
                    if not pending:
                        state['status'] = 'waiting_for_memory'
                        save(root / 'state.json', state)
                        time.sleep(30)
                    break
                row = next(iterator, None)
                if row is None:
                    exhausted = True
                    break
                pending[pool.submit(run_task, row, plan, plan_hash)] = row['membership_sha256']
            if not pending:
                continue
            done, _ = futures.wait(pending, timeout=30, return_when=futures.FIRST_COMPLETED)
            for future in done:
                key = pending.pop(future)
                try:
                    future.result()
                    state['completed'] += 1
                except Exception as error:
                    state['failures'].append(dict(family=key, error=repr(error)))
            state.update(status='running', active_families=list(pending.values()), updated_at=time.time())
            save(root / 'state.json', state)
    state['status'] = 'incomplete_tasks_require_review' if state['failures'] else 'complete_family_jobs_requires_independent_readback'
    save(root / 'state.json', state)


if __name__ == '__main__':
    main()
