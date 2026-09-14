#!/usr/bin/env python3
"""Retain exact baseline ASA rows for unchanged sites in FCS omission inputs."""
import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['baseline', 'baseline-inputs', 'inputs', 'output']:
        p.add_argument('--' + key, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    br = checked_receipt(a.baseline)
    ir = checked_receipt(a.inputs)
    checked_receipt(a.baseline_inputs)
    if br['status'] != 'complete_audited_accessibility_projection':
        raise ValueError('Baseline projection incomplete')
    if br['source_receipts']['inputs'] != sha(a.baseline_inputs / 'receipt.json') or ir['source_input_receipt_sha256'] != sha(a.baseline_inputs / 'receipt.json'):
        raise ValueError('Baseline/sensitivity linkage differs')
    ready = [r for r in read_table(a.inputs / 'marker_summary.tsv') if r['status'] == 'ready_for_inference']
    if len(ready) != ir['ready_markers']:
        raise ValueError('Marker grid differs')
    sequences, columns, seen = {}, {}, {}
    summaries = {}
    for row in ready:
        marker = row['marker']
        folder = a.inputs / marker
        aa = {r.id: str(r.seq) for r in SeqIO.parse(folder / 'aa.faa', 'fasta')}
        di = {r.id: str(r.seq) for r in SeqIO.parse(folder / '3di.faa', 'fasta')}
        base = {r.id: str(r.seq) for r in SeqIO.parse(a.baseline_inputs / marker / 'aa.faa', 'fasta')}
        base_di = {r.id: str(r.seq) for r in SeqIO.parse(a.baseline_inputs / marker / '3di.faa', 'fasta')}
        cs = read_table(folder / 'columns.tsv')
        original_cs = read_table(a.baseline_inputs / marker / 'columns.tsv')
        if set(aa) != set(di) or len(cs) != len(original_cs):
            raise ValueError('Cannot reuse rows if columns changed')
        for before, after in zip(original_cs, cs):
            if any(before[k] != after[k] for k in before) or after['baseline_paired_column_1based'] != after['paired_column_1based']:
                raise ValueError('Column mapping changed')
        columns[marker] = [r['matrix_column_1based'] for r in cs]
        for taxon in aa:
            if aa[taxon] != base[taxon] or di[taxon] != base_di[taxon]:
                raise ValueError('Retained sequence or mask changed')
            sequences[marker, taxon] = aa[taxon], di[taxon]
            seen[marker, taxon] = bytearray(len(cs))
        observed = sum(sum(c != '?' for c in s) for s in aa.values())
        summaries[marker] = {'marker': marker, 'taxa': len(aa), 'paired_columns': len(cs),
                             'observed_sites_linked': observed,
                             'masked_cells_omitted': len(aa) * len(cs) - observed}
    a.output.mkdir(parents=True)
    counts = Counter()
    models = set()
    source_rows = 0
    with gzip.open(a.baseline / 'paired_site_accessibility.tsv.gz', 'rt') as src, gzip.open(a.output / 'paired_site_accessibility.tsv.gz', 'wt') as dst:
        reader = csv.DictReader(src, delimiter='\t')
        writer = csv.DictWriter(dst, reader.fieldnames, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in reader:
            source_rows += 1
            key = row['marker'], row['taxon_id']
            if key not in sequences:
                continue
            index = int(row['paired_column_1based']) - 1
            aa, di = sequences[key]
            if not 0 <= index < len(aa) or seen[key][index] or aa[index] == '?':
                raise ValueError('Duplicate, missing or invalid source cell')
            if aa[index] != row['amino_acid'] or di[index] != row['3di_state'] or columns[key[0]][index] != row['matrix_column_1based']:
                raise ValueError('Source row does not match sensitivity alignment')
            seen[key][index] = 1
            writer.writerow(row)
            counts[key[0]] += 1
            models.add(row['model_id'])
    if source_rows != br['observed_sites_linked']:
        raise ValueError('Baseline count changed')
    for key, bitmap in seen.items():
        if list(bitmap) != [int(c != '?') for c in sequences[key][0]]:
            raise ValueError('Incomplete retained cell grid')
    if any(counts[m] != r['observed_sites_linked'] for m, r in summaries.items()):
        raise ValueError('Counts differ')
    write_table(a.output / 'marker_summary.tsv', list(summaries.values()))
    source_receipts = dict(br['source_receipts'], inputs=sha(a.inputs / 'receipt.json'))
    result = {'status': 'complete_audited_accessibility_projection', 'markers': len(ready),
              'models_used': len(models), 'observed_sites_linked': sum(counts.values()),
              'source_receipts': source_receipts, 'baseline_projection_receipt_sha256': sha(a.baseline / 'receipt.json'),
              'baseline_input_receipt_sha256': sha(a.baseline_inputs / 'receipt.json'),
              'script_sha256': sha(Path(__file__)),
              'interpretation': 'Exact baseline ASA rows retained for every observed FCS sensitivity alignment cell. No coordinate or ASA recomputation; unchanged columns and sequences verified. Full source-residue readback still required for this derived projection.',
              'artifacts': {p.name: sha(p) for p in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
