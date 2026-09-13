#!/usr/bin/env python3
"""Extract native coordinate-derived Foldseek 3Di states and feature descriptors."""
import argparse
import json
import shutil
import subprocess
import urllib.request
from pathlib import Path
from compare_marker_structures import ROOT, sha
from assess_pae_sensitivity import checked_receipt

VERSION = 'e3fadcd07f971e864c094ac4f3a78bf4ed845e07'
SOURCES = ['lib/3di/structureto3di.cpp', 'lib/3di/structureto3di.h',
           'src/strucclustutils/structureto3didescriptor.cpp', 'src/strucclustutils/structcreatedb.cpp']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checked_receipt(args.snapshot)
    if args.output.exists():
        raise FileExistsError('Use a new native extraction output directory')
    executable = Path(shutil.which('foldseek'))
    version = subprocess.check_output([str(executable), 'version'], text=True).strip()
    if version != VERSION:
        raise ValueError('Different Foldseek build requires separate source/validity review')
    models = json.loads((args.snapshot / 'model_provenance.json').read_text())
    names = set()
    for row in models:
        path = ROOT / row['path']
        if sha(path) != row['sha256'] or path.name in names:
            raise ValueError('Changed or ambiguously named model')
        names.add(path.name)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True)
    inputs = args.output / 'inputs'
    inputs.mkdir()
    for row in models:
        path = ROOT / row['path']
        (inputs / path.name).symlink_to(path)
    (args.output / 'model_provenance.json').write_text(json.dumps(models, indent=2) + '\n')
    source = args.output / 'source'
    source.mkdir()
    for name in SOURCES:
        with urllib.request.urlopen(f'https://raw.githubusercontent.com/steineggerlab/foldseek/{version}/{name}', timeout=90) as response:
            (source / Path(name).name).write_bytes(response.read())
    p = args.output
    commands = [[str(executable), 'createdb', str(inputs), str(p / 'structures'), '--threads', '4',
                 '--coord-store-mode', '1', '--mask-bfactor-threshold', '0'],
                [str(executable), 'convert2fasta', str(p / 'structures'), str(p / 'amino_acids.faa')],
                [str(executable), 'lndb', str(p / 'structures_h'), str(p / 'structures_ss_h')],
                [str(executable), 'convert2fasta', str(p / 'structures_ss'), str(p / 'states_3di.faa')],
                [str(executable), 'structureto3didescriptor', str(inputs), str(p / 'descriptors.tsv'), '--threads', '4']]
    config = {'version': version, 'executable_sha256': sha(executable), 'commands': commands,
        'mapping_receipt_sha256': sha(args.snapshot / 'receipt.json'), 'models': len(models),
        'script_sha256': sha(Path(__file__)), 'source_sha256': {p.name: sha(p) for p in source.iterdir()}}
    (p / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
    for i, command in enumerate(commands):
        with (p / f'command-{i}.log').open('w') as handle:
            subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
    print('Native extraction complete; run audit_3di_features.py before using states')


if __name__ == '__main__':
    main()
