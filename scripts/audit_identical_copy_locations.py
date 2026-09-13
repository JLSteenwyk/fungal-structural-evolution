#!/usr/bin/env python3
"""Locate every identical-protein pair assigned to distinct annotated genes."""
import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from audit_busco_gene_copies import ROOT, sha, read_table
from assess_pae_sensitivity import checked_receipt
from map_proteins_to_genes import attributes


def locations(path, mode, wanted):
    spans = defaultdict(set)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('##FASTA'): break
            if not line.strip() or line.startswith('#'): continue
            f = line.rstrip('\n').split('\t')
            if mode == 'creolimax_gtf':
                if f[2] != 'CDS': continue
                match = re.search(r'(?:^|;)\s*gene_id "([^"]+)"', f[8])
                if not match: raise ValueError('Missing exact CDS gene_id')
                gene = match.group(1)
            else:
                if f[2] not in ['gene', 'pseudogene']: continue
                ids = attributes(f[8]).get('ID', [])
                if len(ids) != 1: raise ValueError('Non-unique gene feature ID')
                gene = ids[0]
            if gene in wanted:
                start, end = int(f[3]), int(f[4])
                if not 1 <= start <= end: raise ValueError('Invalid annotation interval')
                spans[gene].add((f[0], start, end, f[6]))
    if mode == 'creolimax_gtf':
        for gene, parts in list(spans.items()):
            if len({(x[0], x[3]) for x in parts}) == 1:
                sample = next(iter(parts))
                spans[gene] = {(sample[0], min(x[1] for x in parts), max(x[2] for x in parts), sample[3])}
    return spans


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sequences', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): raise FileExistsError('Use a new immutable location audit')
    parent = checked_receipt(args.sequences)
    if parent['status'] != 'complete_duplicate_marker_exact_sequence_audit': raise ValueError('Completed sequence audit required')
    pairs = [r for r in read_table(args.sequences / 'protein_pairs.tsv') if r['gene_relationship'] == 'distinct_annotated_genes' and r['exact_full_sequence_identity'] == 'True']
    mapping_path = ROOT / 'metadata/gene_mapping_snapshot.json'
    mapping = {r['taxon_id']: r for r in json.loads(mapping_path.read_text())['taxa']}
    rep_path = ROOT / 'results/gene_representatives/full-v2/receipt.json'
    rep = json.loads(rep_path.read_text())
    if rep['mapping_snapshot_sha256'] != sha(mapping_path): raise ValueError('Gene mapping snapshot differs from representative source')
    reps = {r['taxon_id']: r for r in rep['taxa']}
    ncbi_path = ROOT / 'metadata/annotation_download_receipts.json'
    external_path = ROOT / 'metadata/external_genome_receipts.json'
    annotations = json.loads(ncbi_path.read_text()) + json.loads(external_path.read_text())
    grouped = defaultdict(list)
    for row in pairs: grouped[row['taxon_id']].append(row)
    results, sources = [], []
    for taxon, rows in sorted(grouped.items()):
        info = mapping[taxon]
        if info['mapping_sha256'] != reps[taxon]['mapping_sha256']: raise ValueError('Gene mapping identity differs')
        matches = [r for r in annotations if r.get('sha256') == info['annotation_sha256']]
        if len({r['path'] for r in matches}) != 1: raise ValueError('Cannot resolve exact annotation source')
        path = ROOT / matches[0]['path']
        if sha(path) != info['annotation_sha256']: raise ValueError('Changed annotation')
        wanted = {g for r in rows for key in ['gene_ids_a_json', 'gene_ids_b_json'] for g in json.loads(r[key])}
        spans = locations(path, info['mapping_mode'], wanted)
        for row in rows:
            ga, gb = json.loads(row['gene_ids_a_json']), json.loads(row['gene_ids_b_json'])
            if len(ga) != 1 or len(gb) != 1 or ga == gb: raise ValueError('Distinct unique genes required')
            a, b = sorted(spans.get(ga[0], [])), sorted(spans.get(gb[0], []))
            overlap, gap, strand = '', '', ''
            if len(a) != 1 or len(b) != 1:
                relation = 'missing_or_ambiguous_gene_interval'
            elif a[0][0] != b[0][0]:
                relation = 'different_assembly_sequences'
            else:
                x, y = a[0], b[0]
                overlap = max(0, min(x[2], y[2]) - max(x[1], y[1]) + 1)
                gap = 0 if overlap else max(x[1], y[1]) - min(x[2], y[2]) - 1
                strand = x[3] == y[3] if x[3] in ['+', '-'] and y[3] in ['+', '-'] else ''
                relation = 'same_sequence_overlapping_intervals' if overlap else 'same_sequence_disjoint_intervals'
            results.append(dict(row, annotation_mode=info['mapping_mode'],
                interval_basis='CDS_bounding_span' if info['mapping_mode'] == 'creolimax_gtf' else 'annotated_gene_feature',
                intervals_a_json=json.dumps(a), intervals_b_json=json.dumps(b), location_relationship=relation,
                overlap_bases=overlap, intervening_bases=gap, same_strand=strand))
        sources.append({'taxon_id': taxon, 'annotation_path': str(path.relative_to(ROOT)), 'annotation_sha256': sha(path)})
    if len(results) != parent['exact_identical_distinct_gene_pairs']: raise ValueError('Pair coverage differs')
    args.output.mkdir(parents=True)
    table = args.output / 'identical_copy_locations.tsv'
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, list(results[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(results)
    result = {'status': 'complete_identical_copy_annotation_location_audit', 'pairs': len(results), 'taxa': len(grouped),
        'location_counts': dict(Counter(r['location_relationship'] for r in results)),
        'interval_basis_counts': dict(Counter(r['interval_basis'] for r in results)),
        'sequence_audit_receipt_sha256': sha(args.sequences / 'receipt.json'), 'mapping_snapshot_sha256': sha(mapping_path),
        'representative_receipt_sha256': sha(rep_path), 'script_sha256': sha(Path(__file__)), 'sources': sources,
        'interpretation': 'Coordinates are 1-based inclusive annotation intervals. Assembly sequence identifiers need not represent separate chromosomes. Creolimax uses exact CDS gene_id bounding spans because its gene-feature names differ. Different assembly sequences do not prove haplotigs, and disjoint or overlapping gene intervals do not alone establish biological duplication or annotation error. No copies are removed.',
        'artifacts': {table.name: sha(table)}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'sources'}, indent=2))


if __name__ == '__main__': main()
