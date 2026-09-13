#!/usr/bin/env python3
"""Search all pinned Pfam profiles against the full marker protein dataset."""
import argparse
import concurrent.futures
import fcntl
import json
import shutil
import subprocess
import time
from pathlib import Path
from prepare_pfam import ROOT, digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    pfam_path = ROOT / 'metadata/pfam_release_receipt.json'
    pfam = json.loads(pfam_path.read_text())
    inputs = ROOT / 'data/domains/marker-inputs-v1'
    input_receipt = json.loads((inputs / 'receipt.json').read_text())
    for name, expected in input_receipt['artifacts'].items():
        if digest(inputs / name) != expected:
            raise ValueError('Changed domain search input')
    if sum(r['profiles'] for r in pfam['chunks']) != pfam['families']:
        raise ValueError('Incomplete Pfam profile partition')
    for chunk in pfam['chunks']:
        if digest(ROOT / chunk['path']) != chunk['sha256']:
            raise ValueError('Changed Pfam profile chunk')
    executable = shutil.which('hmmsearch')
    if executable is None:
        raise FileNotFoundError('hmmsearch')
    version = subprocess.run([executable, '-h'], capture_output=True, text=True, check=True).stdout.splitlines()[1]
    if 'HMMER 3.4' not in version:
        raise ValueError('Expected HMMER 3.4')
    args.output.mkdir(parents=True, exist_ok=True)
    lock = (args.output / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    config = {'pfam_receipt_sha256': digest(pfam_path), 'input_receipt_sha256': digest(inputs / 'receipt.json'),
              'executable': executable, 'executable_sha256': digest(Path(executable)), 'version': version,
              'workers': 4, 'hmmer_worker_threads_per_job': 2, 'thresholds': 'Pfam sequence and domain gathering thresholds (--cut_ga)',
              'search_direction': 'hmmsearch: target=protein, query=Pfam model; E-values use unique marker-sequence database size'}
    config_path = args.output / 'config.json'
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        raise ValueError('Changed search configuration; use a new output directory')
    config_path.write_text(json.dumps(config, indent=2) + '\n')

    def search(chunk):
        key = Path(chunk['path']).stem
        receipt_path = args.output / (key + '.receipt.json')
        table = args.output / (key + '.domtblout')
        command = [executable, '--cut_ga', '--cpu', '2', '--noali', '--domtblout', str(table),
                   '-o', '/dev/null', str(ROOT / chunk['path']), str(inputs / 'sequences.faa')]
        if receipt_path.exists():
            r = json.loads(receipt_path.read_text())
            if (r['command'] != command or r['config_sha256'] != digest(config_path)
                    or r['profile_sha256'] != chunk['sha256'] or digest(table) != r['table_sha256']):
                raise ValueError('Stale completed domain search')
            return r
        if table.exists():
            raise FileExistsError('Unreceipted search output requires review before restart: ' + str(table))
        start = time.monotonic()
        with (args.output / (key + '.stderr.log')).open('w') as log:
            subprocess.run(command, stdout=log, stderr=log, check=True)
        tail = table.read_bytes()[-100:]
        if b'# [ok]' not in tail:
            raise ValueError('HMMER domain table lacks completion marker')
        hits = 0
        with table.open() as handle:
            for line in handle:
                if line.strip() and not line.startswith('#'):
                    fields = line.split(maxsplit=22)
                    if len(fields) < 22 or not fields[0].startswith('S') or not fields[4].startswith('PF'):
                        raise ValueError('Unexpected hmmsearch domain table orientation')
                    hits += 1
        r = {'chunk': key, 'profiles': chunk['profiles'], 'domain_hit_rows': hits,
             'elapsed_seconds': time.monotonic() - start, 'command': command,
             'config_sha256': digest(config_path), 'profile_sha256': chunk['sha256'], 'table_sha256': digest(table)}
        receipt_path.write_text(json.dumps(r, indent=2) + '\n')
        return r

    records = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(search, chunk) for chunk in pfam['chunks']]
        try:
            for future in concurrent.futures.as_completed(futures):
                r = future.result()
                records.append(r)
                print(r['chunk'], r['domain_hit_rows'], round(r['elapsed_seconds'], 1), flush=True)
        except Exception:
            for future in futures:
                future.cancel()
            raise
    result = {'status': 'complete_raw_domain_search', 'pfam_families': pfam['families'],
              'input_unique_sequences': input_receipt['unique_sequences'], 'input_marker_proteins': input_receipt['marker_proteins'],
              'domain_hit_rows': sum(r['domain_hit_rows'] for r in records), 'chunks': sorted(records, key=lambda r: r['chunk']),
              'config_sha256': digest(config_path), 'script_sha256': digest(Path(__file__)),
              'interpretation': 'Raw gathering-threshold matches; overlapping hits, Pfam types and domain architectures need interpretation. No hit is not proof of domain absence.'}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print('Completed all', pfam['families'], 'Pfam profiles', flush=True)


if __name__ == '__main__':
    main()
