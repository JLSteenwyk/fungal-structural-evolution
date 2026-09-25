#!/usr/bin/env python3
"""Locate noncanonical residues in checksum-bound residual marker sequences."""
import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from pathlib import Path

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--gap-plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sources = {}

    def read(path, expected=None):
        path = Path(path)
        raw = path.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if expected is not None and digest != expected:
            raise ValueError('Changed input: ' + str(path))
        sources[str(path)] = digest
        return raw

    summary = json.loads(read(args.summary))
    plan = json.loads(read(args.gap_plan))
    tables = {}
    for path, digest in summary['source_sha256'].items():
        raw = read(path, digest)
        if path.endswith('.tsv'):
            tables[Path(path).name] = list(csv.DictReader(io.StringIO(raw.decode()), delimiter='\t'))
    residual = {r['sequence_sha256'] for r in tables['gap_sequence_availability.tsv']
                if int(r['verified_cached_candidates']) == 0}
    targets = {r['sequence_sha256']: r for r in tables['gap_unique_sequences.tsv']
               if r['sequence_sha256'] in residual and r['noncanonical_symbols']}
    links = [r for r in tables['gap_marker_records.tsv'] if r['sequence_sha256'] in targets]
    observed = {}
    canonical = set('ACDEFGHIKLMNPQRSTVWY')
    for marker in sorted({r['marker'] for r in links}):
        path = str(Path(plan['unaligned']) / (marker + '.faa'))
        raw = read(path, plan['pins'][path])
        for record in SeqIO.parse(io.StringIO(raw.decode()), 'fasta'):
            sequence = str(record.seq)
            digest = hashlib.sha256(sequence.encode()).hexdigest()
            if digest not in targets:
                continue
            positions = [[i, aa] for i, aa in enumerate(sequence, 1) if aa not in canonical]
            row = targets[digest]
            if len(sequence) != int(row['length']) or ''.join(sorted(set(a for _, a in positions))) != row['noncanonical_symbols']:
                raise ValueError('Sequence characterization mismatch: ' + digest)
            observed[digest] = dict(length=len(sequence), positions_1based=positions,
                                    symbol_counts=dict(Counter(a for _, a in positions)))
    if set(observed) != set(targets):
        raise ValueError('Missing source sequences')
    symbol_sequences = Counter()
    symbol_residues = Counter()
    for row in observed.values():
        symbol_sequences.update(row['symbol_counts'].keys())
        symbol_residues.update(row['symbol_counts'])
    result = dict(status='complete_residual_marker_symbol_inspection',
                  unique_sequences=len(observed), marker_records=len(links),
                  sequences_per_symbol=dict(symbol_sequences), residues_per_symbol=dict(symbol_residues),
                  sequences=observed, links=links, source_sha256=sources,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope='Exact original sequence symbols and 1-based positions. No sequence correction, imputation, prediction, or inference of the cause of ambiguity.')
    for path, digest in list(sources.items()):
        read(path, digest)
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'unique_sequences', 'marker_records', 'sequences_per_symbol', 'residues_per_symbol']}, indent=2))


if __name__ == '__main__':
    main()
