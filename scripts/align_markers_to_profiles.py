#!/usr/bin/env python3
"""Align the full marker panel to pinned BUSCO profiles; retain match-state columns."""
import argparse
import csv
import fcntl
import hashlib
import json
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from Bio import AlignIO, SeqIO

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--markers', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    extraction = json.loads((args.markers / 'receipt.json').read_text())
    if extraction['status'] != 'complete_extraction':
        raise ValueError('Full marker extraction required')
    executable = shutil.which('hmmalign')
    if executable is None:
        raise FileNotFoundError('hmmalign')
    version = subprocess.run([executable, '-h'], capture_output=True, text=True, check=True).stdout.splitlines()[1]
    with (args.markers / 'occupancy.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    args.output.mkdir(parents=True, exist_ok=True)
    lock = (args.output / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def run(row):
        marker = row['marker']
        source = (args.markers / 'unaligned' / f'{marker}.faa').resolve()
        profile = ROOT / 'data/busco_downloads/lineages/eukaryota_odb12.2/hmms' / f'{marker}.hmm'
        if sha(source) != row['sha256']:
            raise ValueError('Changed source marker')
        stockholm = args.output / f'{marker}.sto'
        target = args.output / f'{marker}.faa'
        receipt_path = args.output / f'{marker}.receipt.json'
        command = [executable, '--amino', '-o', str(stockholm), str(profile), str(source)]
        if receipt_path.exists():
            old = json.loads(receipt_path.read_text())
            if (old['input_sha256'] == sha(source) and old['profile_sha256'] == sha(profile)
                    and old['output_sha256'] == sha(target) and old['stockholm_sha256'] == sha(stockholm)
                    and old['version'] == version):
                return old
            raise ValueError('Stale profile alignment; use a new output directory')
        if target.exists() or stockholm.exists():
            raise FileExistsError('Unreceipted profile output requires review')
        with (args.output / f'{marker}.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        alignment = AlignIO.read(stockholm, 'stockholm')
        reference = alignment.column_annotations['reference_annotation']
        columns = [i for i, symbol in enumerate(reference) if symbol not in '.- ']
        expected_length = int(next(line.split()[1] for line in profile.read_text().splitlines() if line.startswith('LENG ')))
        if len(columns) != expected_length:
            raise ValueError('Profile match-state count differs from HMM length')
        with source.open() as handle:
            originals = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        if len(alignment) != len(originals) or {r.id for r in alignment} != set(originals):
            raise ValueError('Profile alignment changed taxon identities')
        with target.open('w') as out:
            for record in alignment:
                seq = str(record.seq).upper().replace('.', '-')
                if seq.replace('-', '') != str(originals[record.id].seq).upper():
                    raise ValueError('Full Stockholm alignment changed input residues')
                out.write(f'>{record.id}\n' + ''.join(seq[i] for i in columns) + '\n')
        result = {'marker': marker, 'status': 'aligned', 'method': 'HMMER_profile_match_states',
                  'input_sha256': sha(source), 'profile_sha256': sha(profile),
                  'output_sha256': sha(target), 'stockholm_sha256': sha(stockholm),
                  'command': command, 'version': version, 'taxa': len(alignment), 'sites': len(columns),
                  'retained_stockholm_columns_1based': [i + 1 for i in columns]}
        receipt_path.write_text(json.dumps(result, indent=2) + '\n')
        return result

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(run, row) for row in rows]):
            result = future.result()
            results.append(result)
            print(result['marker'], result['sites'], flush=True)
    (args.output / 'receipt.json').write_text(json.dumps({
        'marker_receipt_sha256': sha(args.markers / 'receipt.json'),
        'alignments': sorted(results, key=lambda r: r['marker'])}, indent=2) + '\n')


if __name__ == '__main__':
    main()
