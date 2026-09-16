#!/usr/bin/env python3
"""Describe complete family sequence composition, redundancy and alignment coverage."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
import pandas as pd
from Bio import SeqIO


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    a = ap.parse_args(); plan = json.loads(a.plan.read_text())
    for p, expected in plan['pins'].items():
        if sha(Path(p)) != expected:
            raise ValueError('Changed pinned source ' + p)
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    start = time.time()
    with Path(plan['manifest']).open() as handle:
        manifest = {x['taxon_id']: x for x in csv.DictReader(handle, delimiter='\t')}
    mapping = {}
    for line in Path(plan['species_ids']).read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            key, value = line.split(': ', 1); mapping[key] = value.rsplit('.', 1)[0]
    canonical = set('ACDEFGHIKLMNPQRSTVWY')
    lookup = np.zeros(256, dtype=bool)
    lookup[[ord(c) for c in canonical]] = True
    source = {}; records = {}; global_hashes = Counter(); taxon_hashes = Counter()
    entropy_error = 0.
    for record in SeqIO.parse(plan['source'], 'fasta'):
        sequence = str(record.seq)
        if not sequence or record.id in source:
            raise ValueError('Empty or duplicate source identifier')
        native = record.id.split('_', 1)[0]; taxon = mapping[native]
        if taxon not in manifest:
            raise ValueError('Unknown source taxon')
        counts = Counter(c for c in sequence if c in canonical); n = sum(counts.values())
        if not n:
            raise ValueError('No canonical residues in source sequence')
        entropy = -sum((v / n) * math.log2(v / n) for v in counts.values())
        alternate = math.log2(n) - sum(v * math.log2(v) for v in counts.values()) / n
        entropy_error = max(entropy_error, abs(entropy - alternate))
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        global_hashes[digest] += 1; taxon_hashes[taxon, digest] += 1
        source[record.id] = sequence
        records[record.id] = dict(gene_id=record.id, native_species_id=native, taxon_id=taxon,
                                 sequence_sha256=digest, sequence_length=len(sequence),
                                 canonical_fraction=n / len(sequence), canonical_entropy_bits=entropy,
                                 effective_residue_alphabet=2**entropy,
                                 maximum_canonical_residue_fraction=max(counts.values()) / n)
    if len(records) != plan['expected_genes'] or entropy_error > 1e-12:
        raise ValueError('Source scope or entropy arithmetic differs')
    focal = plan['focal_taxon']
    counts = Counter(x['taxon_id'] for x in records.values())
    if counts[focal] != plan['expected_focal_genes']:
        raise ValueError('Focal taxon count differs')
    columns = plan['expected_alignment_columns']
    occupancy = np.zeros((2, columns), dtype=np.int64)
    canonical_occupancy = np.zeros_like(occupancy)
    seen = set()
    for alignment in SeqIO.parse(plan['alignment'], 'fasta'):
        gene = alignment.id; text = str(alignment.seq)
        if gene in seen or gene not in records or len(text) != columns:
            raise ValueError('Alignment identity or dimension mismatch')
        seen.add(gene)
        retained = text.replace('-', '')
        # Independently recheck ordered residue preservation for this diagnostic.
        original = iter(source[gene])
        if not all(any(c == residue for c in original) for residue in retained):
            raise ValueError('Aligned residues do not preserve source order')
        arr = np.frombuffer(text.encode('ascii'), dtype=np.uint8)
        group = 0 if records[gene]['taxon_id'] == focal else 1
        occupancy[group] += arr != ord('-')
        canonical_occupancy[group] += lookup[arr]
        records[gene].update(alignment_retained_residues=len(retained),
                             source_retention_fraction=len(retained) / len(source[gene]),
                             alignment_nongap_fraction=len(retained) / columns,
                             global_identical_sequence_count=global_hashes[records[gene]['sequence_sha256']],
                             within_taxon_identical_sequence_count=taxon_hashes[records[gene]['taxon_id'], records[gene]['sequence_sha256']])
    if seen != set(records):
        raise ValueError('Alignment/source universe differs')
    frame = pd.DataFrame(records.values()).sort_values('gene_id')
    frame.to_csv(out / 'protein_quality.tsv', sep='\t', index=False)
    sites = pd.DataFrame({'alignment_column_1based': np.arange(1, columns + 1),
                          'focal_nongap': occupancy[0], 'other_nongap': occupancy[1],
                          'focal_canonical': canonical_occupancy[0], 'other_canonical': canonical_occupancy[1]})
    sites.to_csv(out / 'alignment_column_coverage.tsv', sep='\t', index=False)
    assert int(occupancy.sum()) == int(frame.alignment_retained_residues.sum())
    summaries = []
    for label, subset in [('all', frame), ('focal', frame[frame.taxon_id == focal]), ('other', frame[frame.taxon_id != focal])]:
        row = dict(group=label, genes=len(subset), taxa=subset.taxon_id.nunique(),
                   distinct_sequences=subset.sequence_sha256.nunique(),
                   genes_with_within_taxon_exact_duplicate=int((subset.within_taxon_identical_sequence_count > 1).sum()))
        for metric in ['sequence_length', 'canonical_fraction', 'canonical_entropy_bits', 'effective_residue_alphabet', 'maximum_canonical_residue_fraction', 'source_retention_fraction', 'alignment_nongap_fraction']:
            for q, value in subset[metric].quantile([0, .1, .5, .9, 1]).items():
                row[f'{metric}_q{q:g}'] = float(value)
        summaries.append(row)
    pd.DataFrame(summaries).to_csv(out / 'group_summary.tsv', sep='\t', index=False)
    for p, expected in plan['pins'].items():
        if sha(Path(p)) != expected:
            raise ValueError('Source changed during profiling')
    receipt = dict(status='complete_large_family_sequence_quality_profile', family=plan['family'],
                   genes=len(frame), taxa=len(counts), focal_taxon=focal, focal_genes=counts[focal],
                   alignment_columns=columns, entropy_arithmetic_max_difference=entropy_error,
                   retained_residues=int(occupancy.sum()), distinct_sequences=len(global_hashes),
                   elapsed_seconds=time.time() - start, plan_sha256=sha(a.plan), script_sha256=sha(Path(__file__)),
                   artifacts={p.name: sha(p) for p in out.iterdir()},
                   interpretation='All family members retained. Composition quantiles, exact-sequence redundancy and non-gap/canonical alignment coverage are descriptive annotation/homology diagnostics. No low-complexity threshold, gene exclusion, transposon assignment, functional classification or biological duplication is inferred. Original inputs and live tree repair are unchanged.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(pd.DataFrame(summaries)[['group', 'genes', 'distinct_sequences', 'sequence_length_q0.5', 'canonical_entropy_bits_q0.5', 'alignment_nongap_fraction_q0.5', 'source_retention_fraction_q0.5']].to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
