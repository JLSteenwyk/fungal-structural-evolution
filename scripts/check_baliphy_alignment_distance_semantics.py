#!/usr/bin/env python3
"""Bind AxA pairwise semantics to residue-coordinate features and exact binary."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random
import subprocess
import tempfile
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def features(alignment):
    """Canonical/gap fixtures only; do not silently assign semantics to '?' or X."""
    if len({len(s) for s in alignment.values()}) != 1:
        raise ValueError('Unequal aligned lengths')
    if any(set(s) - set('ACDEFGHIKLMNPQRSTVWY-') for s in alignment.values()):
        raise ValueError('Only canonical amino acids and gaps are qualified here')
    tips = sorted(alignment)
    positions = dict.fromkeys(tips, 0)
    result = {}
    for column in zip(*(alignment[t] for t in tips)):
        residues = {}
        for tip, letter in zip(tips, column):
            if letter != '-':
                positions[tip] += 1
                residues[tip] = positions[tip]
            else:
                residues[tip] = None
        for a, b in itertools.combinations(tips, 2):
            x, y = residues[a], residues[b]
            if x is not None or y is not None:
                result[a, x, b, y] = int(x is not None) + int(y is not None)
    return result


def distance(a, b):
    if {t: s.replace('-', '') for t, s in a.items()} != {t: s.replace('-', '') for t, s in b.items()}:
        raise ValueError('Extant ungapped identities must agree')
    x, y = features(a), features(b)
    return sum(x[k] for k in x.keys() - y.keys()) + sum(y[k] for k in y.keys() - x.keys())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    binary = Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/alignment-distances').resolve()
    source = Path('data/software_audits/baliphy-alignment-distances-20260928-v1')
    rng = random.Random(20260928)
    alignments = [{'a': 'AA', 'b': 'A-'}, {'a': 'AA', 'b': '-A'}]
    assert distance(*alignments) == 6
    # Full matrix requires the same tip set and ungapped sequences throughout.
    samples = []
    sequences = {'a': 'ACAA', 'b': 'AAA', 'c': 'CA', 'd': 'A'}
    for _ in range(24):
        offsets = dict.fromkeys(sequences, 0)
        a = dict.fromkeys(sequences, '')
        while any(offsets[t] < len(s) for t, s in sequences.items()):
            available = [t for t in sequences if offsets[t] < len(sequences[t])]
            chosen = rng.sample(available, rng.randint(1, len(available)))
            for t, s in sequences.items():
                a[t] += s[offsets[t]] if t in chosen else '-'
                offsets[t] += t in chosen
        samples.append(a)
    # Add duplicate, redundant all-gap column and reordered FASTA record fixtures.
    samples.extend([samples[0], {t: '-' + s for t, s in samples[0].items()}, dict(reversed(list(samples[0].items())))])
    groups = [alignments, samples]
    started = time.monotonic()
    evidence = []
    with tempfile.TemporaryDirectory() as temporary:
        for group_id, group in enumerate(groups):
            paths = []
            for i, alignment in enumerate(group):
                path = Path(temporary) / f'{group_id}-{i}.fasta'
                path.write_text(''.join(f'>{t}\n{s}\n' for t, s in alignment.items()))
                paths.append(str(path))
            command = [str(binary), 'AxA', *paths, '--alphabet', 'Amino-Acids', '--distances', 'pairwise']
            completed = subprocess.run(command, text=True, capture_output=True, check=True, timeout=60)
            actual = [[float(v) for v in line.split()] for line in completed.stdout.splitlines() if line.strip()]
            expected = [[distance(a, b) for b in group] for a in group]
            assert actual == expected, (actual, expected)
            evidence.append(dict(alignments=group, expected=expected, actual=actual, stderr=completed.stderr))
    for invalid in [{'a': 'X'}, {'a': '?'}, {'a': 'A', 'b': 'AA'}]:
        try:
            features(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid fixture accepted')
    pins = {str(p): sha(p) for p in [binary, Path(__file__), *sorted(source.glob('*.cc')), source / 'source.json']}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(dict(status='passed_canonical_fixture_metric_semantics', pins=pins,
        source_commit=json.loads((source / 'source.json').read_text())['commit'],
        matrix_entries=sum(len(g)**2 for g in groups), elapsed_seconds=time.monotonic()-started, evidence=evidence,
        interpretation='AxA pairwise is the weighted symmetric difference of residue-pair and residue-gap features; residue-pair weight 2, residue-gap weight 1. It is unnormalized.',
        limitations='Canonical/gap synthetic fixtures only. No posterior mixing or convergence claim. X/? handling and full production extant projection require separate validation.'), indent=2)+'\n')
    print(json.dumps(dict(matrix_entries=sum(len(g)**2 for g in groups), elapsed_seconds=time.monotonic()-started)))


if __name__ == '__main__':
    main()
