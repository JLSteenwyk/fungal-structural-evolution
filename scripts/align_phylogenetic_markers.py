#!/usr/bin/env python3
"""Align a complete marker extraction with bounded parallelism and verified resume."""
import argparse
import csv
import fcntl
import hashlib
import json
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from Bio import SeqIO


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_alignment(source, aligned):
    with source.open() as handle:
        source_records = list(SeqIO.parse(handle, 'fasta'))
    before = {r.id: str(r.seq).upper() for r in source_records}
    if len(source_records) != len(before):
        raise ValueError('Duplicate input taxon identities')
    with aligned.open() as handle:
        records = list(SeqIO.parse(handle, 'fasta'))
    after = {r.id: str(r.seq).upper() for r in records}
    if not before or len(records) != len(after) or before.keys() != after.keys():
        raise ValueError('Alignment changed taxon identities or duplicated rows')
    if len({len(seq) for seq in after.values()}) != 1:
        raise ValueError('Unequal alignment row lengths')
    if any(after[name].replace('-', '') != seq for name, seq in before.items()):
        raise ValueError('Alignment changed input residues')
    return len(records), len(next(iter(after.values())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--markers', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args()
    if args.workers < 1 or args.threads < 1:
        parser.error('Worker and thread counts must be positive')
    receipt = json.loads((args.markers / 'receipt.json').read_text())
    if receipt['status'] != 'complete_extraction' or receipt['pending_taxa']:
        raise ValueError('Incomplete staging snapshots cannot enter alignment')
    executable = shutil.which('mafft')
    if executable is None:
        raise FileNotFoundError('mafft')
    version = subprocess.run([executable, '--version'], capture_output=True, text=True, check=True)
    version = (version.stdout + version.stderr).strip()
    with (args.markers / 'occupancy.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    args.output.mkdir(parents=True, exist_ok=True)
    lock = (args.output / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def run(row):
        marker = row['marker']
        source = (args.markers / 'unaligned' / f'{marker}.faa').resolve()
        if sha(source) != row['sha256']:
            raise ValueError(f'Changed marker input: {marker}')
        if int(row['single_copy_taxa']) < 4:
            return {'marker': marker, 'status': 'insufficient_taxa'}
        target = args.output / f'{marker}.faa'
        record_path = args.output / f'{marker}.receipt.json'
        command = [executable, '--auto', '--thread', str(args.threads), '--inputorder', str(source)]
        if record_path.exists():
            old = json.loads(record_path.read_text())
            if (old['input_sha256'] == row['sha256'] and old['command'] == command
                    and old['version'] == version and target.exists() and old['output_sha256'] == sha(target)):
                validate_alignment(source, target)
                return old
            raise ValueError(f'Stale alignment receipt: {marker}; select a new output directory')
        if target.exists():
            raise FileExistsError(f'Unreceipted alignment requires review: {target}')
        temp = target.with_suffix('.partial')
        with temp.open('w') as out, (args.output / f'{marker}.log').open('w') as log:
            subprocess.run(command, stdout=out, stderr=log, check=True)
        taxa, sites = validate_alignment(source, temp)
        temp.replace(target)
        result = {'marker': marker, 'status': 'aligned', 'input_sha256': row['sha256'],
                  'output_sha256': sha(target), 'command': command, 'version': version,
                  'taxa': taxa, 'sites': sites}
        record_path.write_text(json.dumps(result, indent=2) + '\n')
        return result

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(run, row) for row in rows]):
            result = future.result()
            results.append(result)
            print(result['marker'], result['status'], flush=True)
    (args.output / 'receipt.json').write_text(json.dumps({
        'marker_receipt_sha256': sha(args.markers / 'receipt.json'),
        'alignments': sorted(results, key=lambda r: r['marker'])}, indent=2) + '\n')


if __name__ == '__main__':
    main()
