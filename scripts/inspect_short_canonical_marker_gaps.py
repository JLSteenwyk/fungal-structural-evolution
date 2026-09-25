#!/usr/bin/env python3
"""Identify residual marker gaps within the existing sequence/length limits."""
import argparse
import csv
import hashlib
import io
import json
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

    def table(raw):
        return list(csv.DictReader(io.StringIO(raw.decode()), delimiter='\t'))

    summary = json.loads(read(args.summary))
    plan = json.loads(read(args.gap_plan))
    tables = {}
    for path, digest in summary['source_sha256'].items():
        raw = read(path, digest)
        if path.endswith('.tsv'):
            tables[Path(path).name] = table(raw)
    residual = {r['sequence_sha256'] for r in tables['gap_sequence_availability.tsv']
                if int(r['verified_cached_candidates']) == 0}
    targets = {r['sequence_sha256']: r for r in tables['gap_unique_sequences.tsv']
               if r['sequence_sha256'] in residual and int(r['length']) <= 1024
               and not r['noncanonical_symbols']}
    if len(targets) != summary['disjoint_length_alphabet_categories']['at_most_1024_canonical']:
        raise ValueError('Summary count mismatch')
    manifest = table(read(plan['manifest'], plan['pins'][plan['manifest']]))
    taxa = {r['taxon_id']: r for r in manifest}
    if len(taxa) != len(manifest):
        raise ValueError('Duplicate taxon')
    rows = []
    for link in tables['gap_marker_records.tsv']:
        digest = link['sequence_sha256']
        if digest not in targets:
            continue
        path = str(Path(plan['unaligned']) / (link['marker'] + '.faa'))
        records = [r for r in SeqIO.parse(io.StringIO(read(path, plan['pins'][path]).decode()), 'fasta')
                   if r.id == link['taxon_id']]
        if len(records) != 1:
            raise ValueError('Nonunique marker/taxon sequence')
        seq = str(records[0].seq)
        if (hashlib.sha256(seq.encode()).hexdigest() != digest
                or len(seq) != int(targets[digest]['length'])
                or not set(seq) <= set('ACDEFGHIKLMNPQRSTVWY')):
            raise ValueError('Sequence mismatch')
        taxon = taxa[link['taxon_id']]
        rows.append(dict(marker=link['marker'], taxon_id=link['taxon_id'],
                         protein_id=link['protein_id'], species_name=taxon['species_name'],
                         assembly_accession=taxon['assembly_accession'],
                         length=len(seq), sequence_sha256=digest))
    if {r['sequence_sha256'] for r in rows} != set(targets):
        raise ValueError('Missing marker link')
    for path, digest in list(sources.items()):
        read(path, digest)
    result = dict(status='complete_short_canonical_residual_inventory',
                  unique_sequences=len(targets), marker_records=len(rows), records=rows,
                  source_sha256=sources,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope='Frozen cache gaps with exact canonical sequences at most 1024 residues. '
                        'Not a GPU launch, a feasibility guarantee, or proof of current database absence.')
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    main()
