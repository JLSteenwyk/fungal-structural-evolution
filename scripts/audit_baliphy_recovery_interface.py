#!/usr/bin/env python3
"""Capture installed recovery interfaces without claiming exhaustive absence."""
import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan_path = Path('metadata/baliphy_prior_initialization_plan_20260927.json')
    plan = json.loads(plan_path.read_text())
    binary = Path(plan['binary'])
    assert sha(binary) == plan['pins'][str(binary)]
    modules = binary.parent.parent / 'lib/bali-phy/haskell'
    sources = [modules / 'MCMC.hs', *sorted((modules / 'MCMC').rglob('*.hs'))]
    assert len(sources) > 3 and all(p.is_file() for p in sources)
    args.output.mkdir(parents=True, exist_ok=False)
    pattern = re.compile(r'checkpoint|restart|resume|serializ|restore', re.I)
    hits = []
    pins = {str(plan_path): sha(plan_path), str(binary): sha(binary),
            str(Path(__file__)): sha(__file__)}
    for level in ['advanced', 'expert', 'developer']:
        command = [str(binary), 'help', level]
        result = subprocess.run(command, text=True, capture_output=True, check=True)
        target = args.output / (level + '.txt')
        target.write_text(result.stdout + result.stderr)
        for i, line in enumerate(target.read_text().splitlines(), 1):
            if pattern.search(line):
                hits.append(dict(path=str(target), line=i, text=line))
    for source in sources:
        pins[str(source)] = sha(source)
        for i, line in enumerate(source.read_text().splitlines(), 1):
            if pattern.search(line):
                hits.append(dict(path=str(source), line=i, text=line))
    receipt = dict(status='completed_bounded_installed_interface_search',
                   searched_modules=len(sources), search_pattern=pattern.pattern,
                   matches=hits, pins=pins,
                   artifacts={p.name: sha(p) for p in args.output.iterdir()},
                   scope='Help levels and installed MCMC Haskell modules only. '
                   'No proof that native code or other interfaces lack checkpoint support. '
                   'No restart or interrupted-chain recovery was executed.')
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(modules=len(sources), matches=len(hits))))


if __name__ == '__main__':
    main()
