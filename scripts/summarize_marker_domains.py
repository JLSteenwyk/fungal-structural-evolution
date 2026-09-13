#!/usr/bin/env python3
"""Validate complete Pfam searches and retain annotated hits and overlap ambiguity."""
import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from Bio import SeqIO
from prepare_pfam import ROOT, digest


def read_metadata(path):
    entries, row = {}, {}
    with path.open() as handle:
        for line in handle:
            if line.startswith('#=GF '):
                _, key, value = line.rstrip().split(maxsplit=2)
                row[key] = value.strip()
            elif line.strip() == '//':
                if not {'AC', 'ID', 'TP', 'ML', 'GA'} <= row.keys() or row['AC'] in entries:
                    raise ValueError('Incomplete or duplicate Pfam metadata')
                ga = [float(v.strip()) for v in row['GA'].split(';') if v.strip()]
                if len(ga) != 2 or not all(math.isfinite(x) for x in ga):
                    raise ValueError('Invalid gathering thresholds')
                row.update(sequence_ga=ga[0], domain_ga=ga[1])
                entries[row['AC']] = row
                row = {}
    if row:
        raise ValueError('Truncated Pfam metadata')
    return entries


def parse_hit(line, metadata, lengths):
    f = line.split(maxsplit=22)
    if len(f) < 22 or f[0] not in lengths or f[4] not in metadata:
        raise ValueError('Unknown target or Pfam accession')
    meta = metadata[f[4]]
    tlen, qlen = int(f[2]), int(f[5])
    hf, ht, af, at, ef, et = map(int, f[15:21])
    if (tlen != lengths[f[0]] or qlen != int(meta['ML']) or f[3] != meta['ID']
            or not 1 <= hf <= ht <= qlen or not 1 <= ef <= af <= at <= et <= tlen):
        raise ValueError('Invalid model identity or domain residue coordinates')
    values = [float(f[i]) for i in [6, 7, 8, 11, 12, 13, 14, 21]]
    if not all(math.isfinite(x) for x in values) or any(float(f[i]) < 0 for i in [6, 11, 12]):
        raise ValueError('Invalid score or E-value')
    if not 0 <= float(f[21]) <= 1:
        raise ValueError('Invalid posterior accuracy')
    if float(f[7]) < meta['sequence_ga'] - .051 or float(f[13]) < meta['domain_ga'] - .051:
        raise ValueError('Hit fails gathering thresholds beyond score-rounding tolerance')
    return {'sequence_id': f[0], 'protein_length': tlen, 'pfam_accession': f[4], 'pfam_name': f[3],
        'pfam_type': meta['TP'], 'pfam_clan': meta.get('CL', ''), 'description': meta.get('DE', ''),
        'hmm_length': qlen, 'hmm_start': hf, 'hmm_end': ht, 'alignment_start': af, 'alignment_end': at,
        'envelope_start': ef, 'envelope_end': et, 'hmm_coverage': (ht - hf + 1) / qlen,
        'sequence_score': float(f[7]), 'domain_score': float(f[13]), 'sequence_ga': meta['sequence_ga'],
        'domain_ga': meta['domain_ga'], 'conditional_evalue': float(f[11]), 'independent_evalue': float(f[12]),
        'posterior_accuracy': float(f[21])}


def overlap_size(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--search', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads((args.search / 'receipt.json').read_text())
    if receipt['status'] != 'complete_raw_domain_search':
        raise ValueError('Complete raw domain search required')
    config = json.loads((args.search / 'config.json').read_text())
    if digest(args.search / 'config.json') != receipt['config_sha256']:
        raise ValueError('Changed search configuration')
    pfam_path = ROOT / 'metadata/pfam_release_receipt.json'
    if digest(pfam_path) != config['pfam_receipt_sha256']:
        raise ValueError('Changed Pfam release receipt')
    pfam = json.loads(pfam_path.read_text())
    if sum(r['profiles'] for r in receipt['chunks']) != pfam['families']:
        raise ValueError('Incomplete profile coverage')
    inputs = ROOT / 'data/domains/marker-inputs-v1'
    if digest(inputs / 'receipt.json') != config['input_receipt_sha256']:
        raise ValueError('Changed source receipt')
    ir = json.loads((inputs / 'receipt.json').read_text())
    for name, checksum in ir['artifacts'].items():
        if digest(inputs / name) != checksum:
            raise ValueError('Changed marker inputs')
    metadata_path = ROOT / 'data/pfam/38.2/Pfam-A.hmm.dat'
    item = next(r for r in pfam['files'] if r['uncompressed_file'] == metadata_path.name)
    if digest(metadata_path) != item['uncompressed_sha256']:
        raise ValueError('Changed Pfam metadata')
    metadata = read_metadata(metadata_path)
    if len(metadata) != pfam['families']:
        raise ValueError('Metadata family count mismatch')
    lengths = {r.id: len(r.seq) for r in SeqIO.parse(inputs / 'sequences.faa', 'fasta')}
    hits = []
    for chunk in receipt['chunks']:
        path = args.search / (chunk['chunk'] + '.domtblout')
        if digest(path) != chunk['table_sha256']:
            raise ValueError('Changed search results')
        with path.open() as handle:
            for line in handle:
                if line.strip() and not line.startswith('#'):
                    hits.append(parse_hit(line, metadata, lengths))
    hits.sort(key=lambda r: (r['sequence_id'], r['alignment_start'], r['alignment_end'], r['pfam_accession']))
    by_sequence = defaultdict(list)
    for i, hit in enumerate(hits, 1):
        hit['hit_id'] = f'H{i:08}'
        by_sequence[hit['sequence_id']].append(hit)
    overlaps, ambiguous = [], set()
    for sequence, seq_hits in by_sequence.items():
        for i, a in enumerate(seq_hits):
            for b in seq_hits[i + 1:]:
                overlap = overlap_size(a['alignment_start'], a['alignment_end'], b['alignment_start'], b['alignment_end'])
                if overlap:
                    overlaps.append({'sequence_id': sequence, 'hit_a': a['hit_id'], 'hit_b': b['hit_id'],
                        'alignment_overlap_residues': overlap,
                        'same_clan': bool(a['pfam_clan']) and a['pfam_clan'] == b['pfam_clan'],
                        'same_family': a['pfam_accession'] == b['pfam_accession']})
                    ambiguous.add(sequence)
    summaries = []
    with (inputs / 'protein_links.tsv').open() as handle:
        for link in csv.DictReader(handle, delimiter='\t'):
            rows = by_sequence[link['sequence_id']]
            types = Counter(r['pfam_type'] for r in rows)
            summaries.append({k: link[k] for k in ['marker', 'taxon_id', 'protein_id', 'sequence_id']} | {
                'raw_hit_count': len(rows), 'domain_type_hits': types['Domain'], 'family_type_hits': types['Family'],
                'repeat_type_hits': types['Repeat'], 'overlapping_hits': link['sequence_id'] in ambiguous,
                'status': 'raw_hits_require_architecture_review' if rows else 'no_GA_hit_not_proven_absence'})
    if args.output.exists():
        raise FileExistsError('Use a new immutable annotation snapshot')
    args.output.mkdir(parents=True)
    for filename, rows in [('raw_annotated_hits.tsv', hits), ('overlapping_hits.tsv', overlaps), ('protein_summary.tsv', summaries)]:
        with (args.output / filename).open('w') as handle:
            if rows:
                writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows(rows)
    result = {'search_receipt_sha256': digest(args.search / 'receipt.json'), 'script_sha256': digest(Path(__file__)),
        'unique_sequences': len(lengths), 'unique_sequences_with_hits': len({r['sequence_id'] for r in hits}),
        'marker_proteins': len(summaries), 'raw_hits': len(hits), 'unique_sequences_with_overlapping_hits': len(ambiguous),
        'type_counts': dict(Counter(r['pfam_type'] for r in hits)),
        'interpretation': 'All GA hits retained with overlaps flagged. No resolved architecture, evolutionary gain/loss, orthology or functional validation is implied.',
        'artifacts': {p.name: digest(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
