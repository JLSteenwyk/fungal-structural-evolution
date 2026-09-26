#!/usr/bin/env python3
"""Compare standard, unfiltered and permissive Pfam searches for five focal proteins."""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import time
from pathlib import Path
from Bio import SeqIO
from compare_marker_structures import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable output')
    sequence_root = Path('results/domains/cross-clan-focal-sequences-20260924-v1')
    fasta = sequence_root / 'proteins.faa'
    sequence_receipt = sequence_root / 'receipt.json'
    source = json.loads(sequence_receipt.read_text())
    sequences = list(SeqIO.parse(fasta, 'fasta'))
    expected = {r['label']: r for r in source['proteins']}
    if len(sequences) != 5 or {r.id for r in sequences} != expected.keys():
        raise ValueError('Focal sequence set changed')
    for record in sequences:
        if hashlib.sha256(str(record.seq).encode()).hexdigest() != expected[record.id]['sequence_sha256']:
            raise ValueError('Focal sequence changed')
    pfam_path = Path('metadata/pfam_release_receipt.json')
    pfam = json.loads(pfam_path.read_text())
    entry = next(r for r in pfam['files'] if r['file'] == 'Pfam-A.hmm.gz')
    library = Path('data/pfam/38.2') / entry['uncompressed_file']
    if sha(library) != entry['uncompressed_sha256']:
        raise ValueError('Pfam library changed')
    search, fetch = [shutil.which(name) for name in ['hmmsearch', 'hmmfetch']]
    if not search or not fetch:
        raise FileNotFoundError('HMMER executables required')
    args.output.mkdir(parents=True)
    requested = ['PF00085.27', 'PF26973.1']
    ids = args.output / 'profile_accessions.txt'
    ids.write_text('\n'.join(requested) + '\n')
    profiles = args.output / 'profiles.hmm'
    fetch_command = [fetch, '-f', str(library), str(ids)]
    with profiles.open('w') as handle:
        subprocess.run(fetch_command, stdout=handle, check=True)
    profile_text = profiles.read_text()
    accessions = re.findall(r'^ACC\s+(\S+)', profile_text, re.M)
    if len(accessions) != 2 or set(accessions) != set(requested):
        raise ValueError('Unexpected profile extraction')
    thresholds = {}
    for block in profile_text.split('//'):
        acc = re.search(r'^ACC\s+(\S+)', block, re.M)
        if acc:
            ga = re.search(r'^GA\s+([\d.+-]+)\s+([\d.+-]+);?\s*$', block, re.M)
            if ga is None:
                raise ValueError('Missing gathering thresholds: ' + acc[1])
            thresholds[acc[1]] = {'sequence': float(ga[1]), 'domain': float(ga[2])}
    config = {'profiles': requested, 'gathering_thresholds': thresholds,
              'source_sha256': {str(p): sha(p) for p in [sequence_receipt, fasta, pfam_path, library]},
              'software_sha256': {p: sha(Path(p)) for p in [search, fetch]},
              'fetch_command': fetch_command,
              'version': subprocess.run([search, '-h'], capture_output=True, text=True, check=True).stdout.splitlines()[1],
              'resources': {'worker_threads': 1, 'master_thread_additional': True,
                            'planning_memory_gib': 1, 'planning_wall_minutes': 10,
                            'planning_output_mib': 10, 'gpu': False, 'paid_resources': False},
              'scope': 'Exploratory five-protein/two-profile diagnostic. E-values use five target sequences, not the original proteome search space. Permissive hits are not accepted annotations. No ancestral reconstruction or new domain-event claim.'}
    (args.output / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
    runs = []
    for label, options in [('standard_ga', ['--cut_ga']),
                           ('max_ga', ['--max', '--cut_ga']),
                           ('max_permissive', ['--max', '-T', '-1000', '--domT', '-1000'])]:
        table = args.output / (label + '.domtblout')
        report = args.output / (label + '.txt')
        command = [search, '--cpu', '1', '--seed', '42', '--noali',
                   '--domtblout', str(table), '-o', str(report)] + options + [str(profiles), str(fasta)]
        start = time.monotonic()
        with (args.output / (label + '.stderr.log')).open('w') as handle:
            subprocess.run(command, stdout=handle, stderr=handle, check=True, timeout=600)
        if '[ok]' not in report.read_text()[-100:]:
            raise ValueError('Incomplete HMMER report')
        hits = []
        for line in table.read_text().splitlines():
            if not line or line.startswith('#'):
                continue
            f = line.split(maxsplit=22)
            if len(f) < 22 or f[0] not in expected or f[4] not in requested:
                raise ValueError('Unexpected domain table row')
            ga = thresholds[f[4]]
            hits.append({'protein': f[0], 'profile': f[4], 'sequence_score': float(f[7]),
                         'domain_score': float(f[13]), 'sequence_evalue': float(f[6]),
                         'domain_independent_evalue': float(f[12]),
                         'hmm_start': int(f[15]), 'hmm_end': int(f[16]),
                         'alignment_start': int(f[17]), 'alignment_end': int(f[18]),
                         'envelope_start': int(f[19]), 'envelope_end': int(f[20]),
                         'passes_sequence_and_domain_ga': float(f[7]) >= ga['sequence'] and float(f[13]) >= ga['domain']})
        runs.append({'label': label, 'command': command, 'elapsed_seconds': time.monotonic() - start, 'hits': hits})
    result = {'status': 'complete_focal_profile_search_sensitivity',
              'config_sha256': sha(args.output / 'config.json'), 'script_sha256': sha(Path(__file__)),
              'runs': runs, 'scope': config['scope'],
              'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({r['label']: len(r['hits']) for r in runs}))


if __name__ == '__main__':
    main()
