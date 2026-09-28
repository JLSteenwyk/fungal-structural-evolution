#!/usr/bin/env python3
"""Compile isolated FastML variants from a checksum-bound preparation."""
import argparse
import json
from pathlib import Path
import subprocess
from ancestral_chain_attempt import sha, write_json


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest
    prep = json.loads(Path(plan['preparation']).read_text())
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    builds = []
    for variant in prep['variants']:
        root = Path(variant['directory'])
        for name, digest in prep['source_inventory'].items():
            expected = variant['changes'].get(name, {}).get('modified_sha256', digest)
            assert sha(root / name) == expected, name
        log = out / (variant['label'] + '.log')
        with log.open('w') as handle:
            for command in variant['build_commands']:
                subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
        binary = root / 'programs/gainLoss/gainLoss'
        assert binary.is_file()
        builds.append(dict(label=variant['label'], binary=str(binary), binary_sha256=sha(binary),
                           log=str(log), log_sha256=sha(log)))
    assert all(sha(Path(prep['source']) / name) == digest
               for name, digest in prep['source_inventory'].items())
    write_json(out / 'receipt.json', dict(status='compiled_not_scientifically_qualified', builds=builds,
        plan_sha256=sha(args.plan), preparation_sha256=sha(plan['preparation']),
        compiler=subprocess.check_output(['g++', '--version'], text=True),
        scope='Both isolated variants built successfully; original source and binary unchanged. '
              'Paired inference, independent replay and numerical qualification remain required.'))


if __name__ == '__main__':
    main()
