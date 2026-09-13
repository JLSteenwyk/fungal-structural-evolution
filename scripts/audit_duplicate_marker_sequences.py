#!/usr/bin/env python3
"""Audit exact sequence redundancy within every duplicated raw BUSCO call."""
import argparse
import csv
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path
from Bio import SeqIO
from audit_busco_gene_copies import ROOT, sha, read_table
from assess_pae_sensitivity import checked_receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--copies', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable audit output')
    receipt = checked_receipt(args.copies)
    if receipt['status'] != 'complete_raw_busco_gene_annotation_audit':
        raise ValueError('Complete copy-annotation audit required')
    manifest_path = ROOT / 'metadata/analysis_manifest.tsv'
    if sha(manifest_path) != receipt['manifest_sha256']:
        raise ValueError('Changed sampling manifest')
    manifest = read_table(manifest_path)
    input_path = ROOT / 'metadata/qc_input_receipts.json'
    inputs = {r['taxon_id']: r for r in json.loads(input_path.read_text())}
    reps_path = ROOT / 'results/gene_representatives/full-v2/receipt.json'
    if sha(reps_path) != receipt['representative_receipt_sha256']:
        raise ValueError('Representative provenance changed')
    reps = {r['taxon_id']: r for r in json.loads(reps_path.read_text())['taxa']}
    all_rows = read_table(args.copies / 'marker_gene_copies.tsv')
    by_taxon = defaultdict(list)
    for row in all_rows:
        if row['busco_status'] == 'Duplicated':
            by_taxon[row['taxon_id']].append(row)
    comparisons, markers, sources = [], [], []
    for taxon, calls in sorted(by_taxon.items()):
        source = ROOT / inputs[taxon]['input_path']
        if sha(source) != inputs[taxon]['sha256'] or inputs[taxon]['sha256'] != reps[taxon]['source_sha256']:
            raise ValueError('Changed source protein sequences')
        wanted = {pid for row in calls for pid in json.loads(row['protein_ids_json'])}
        sequences = {}
        for record in SeqIO.parse(source, 'fasta'):
            if record.id in wanted:
                if record.id in sequences:
                    raise ValueError('Repeated source protein identifier')
                sequences[record.id] = str(record.seq)
        if set(sequences) != wanted:
            raise ValueError('Incomplete original BUSCO sequence recovery')
        decision_path = (ROOT / reps[taxon]['path']).with_suffix('.decisions.tsv')
        if sha(decision_path) != reps[taxon]['decisions_sha256']:
            raise ValueError('Changed explicit gene mappings')
        decisions = {r['protein_id']: r for r in read_table(decision_path) if r['protein_id'] in wanted}
        hashes = {pid: hashlib.sha256(seq.encode()).hexdigest() for pid, seq in sequences.items()}
        for call in calls:
            proteins = json.loads(call['protein_ids_json'])
            local = []
            for a, b in itertools.combinations(proteins, 2):
                da, db = decisions[a], decisions[b]
                ga, gb = json.loads(da['gene_ids_json']), json.loads(db['gene_ids_json'])
                if da['status'] == db['status'] == 'unique_gene':
                    if len(ga) != 1 or len(gb) != 1:
                        raise ValueError('Ambiguous unique gene annotation')
                    relation = 'same_annotated_gene' if ga == gb else 'distinct_annotated_genes'
                else:
                    relation = 'unresolved_gene_mapping'
                row = {'taxon_id': taxon, 'marker': call['marker'], 'protein_a': a, 'protein_b': b,
                    'gene_relationship': relation, 'gene_ids_a_json': json.dumps(ga), 'gene_ids_b_json': json.dumps(gb),
                    'length_a': len(sequences[a]), 'length_b': len(sequences[b]),
                    'sequence_sha256_a': hashes[a], 'sequence_sha256_b': hashes[b],
                    'exact_full_sequence_identity': sequences[a] == sequences[b]}
                comparisons.append(row); local.append(row)
            markers.append({'taxon_id': taxon, 'marker': call['marker'], 'hit_proteins': len(proteins),
                'distinct_full_sequences': len({sequences[p] for p in proteins}),
                'duplicate_annotation_category': call['duplicate_annotation_category'],
                'protein_pairs': len(local), 'exact_identical_pairs': sum(r['exact_full_sequence_identity'] for r in local),
                'distinct_gene_pairs': sum(r['gene_relationship'] == 'distinct_annotated_genes' for r in local),
                'exact_identical_distinct_gene_pairs': sum(r['gene_relationship'] == 'distinct_annotated_genes' and r['exact_full_sequence_identity'] for r in local)})
        sources.append({'taxon_id': taxon, 'source_proteome_path': str(source.relative_to(ROOT)),
            'source_proteome_sha256': inputs[taxon]['sha256'], 'decision_sha256': sha(decision_path)})
    grouped = defaultdict(list)
    for row in markers:
        grouped[row['taxon_id']].append(row)
    summaries = []
    for taxon in manifest:
        rows = grouped[taxon['taxon_id']]
        summaries.append({'taxon_id': taxon['taxon_id'], 'species_name': taxon['species_name'],
            'duplicated_markers': len(rows), 'duplicated_markers_one_full_sequence': sum(r['distinct_full_sequences'] == 1 for r in rows),
            'markers_with_identical_distinct_gene_pair': sum(r['exact_identical_distinct_gene_pairs'] > 0 for r in rows),
            'distinct_gene_pairs': sum(r['distinct_gene_pairs'] for r in rows),
            'exact_identical_distinct_gene_pairs': sum(r['exact_identical_distinct_gene_pairs'] for r in rows)})
    if len(markers) != receipt['totals']['raw_duplicated']:
        raise ValueError('Duplicated-call scope differs')
    args.output.mkdir(parents=True)
    for name, rows in [('protein_pairs.tsv', comparisons), ('marker_redundancy.tsv', markers), ('taxon_redundancy.tsv', summaries)]:
        with (args.output / name).open('w') as handle:
            writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
    result = {'status': 'complete_duplicate_marker_exact_sequence_audit', 'taxa_in_scope': len(summaries),
        'taxa_with_duplicated_calls': len(by_taxon), 'duplicated_marker_calls': len(markers), 'protein_pairs': len(comparisons),
        'pair_gene_relationships': dict(Counter(r['gene_relationship'] for r in comparisons)),
        'exact_identical_pairs': sum(r['exact_full_sequence_identity'] for r in comparisons),
        'exact_identical_distinct_gene_pairs': sum(r['exact_identical_distinct_gene_pairs'] for r in markers),
        'markers_with_identical_distinct_gene_pair': sum(r['exact_identical_distinct_gene_pairs'] > 0 for r in markers),
        'copy_audit_receipt_sha256': sha(args.copies / 'receipt.json'), 'qc_input_receipts_sha256': sha(input_path),
        'script_sha256': sha(Path(__file__)), 'sources': sources,
        'interpretation': 'Complete original protein strings compared without alignment, truncation or collapsing loci. Identical proteins from distinct annotated genes are compatible with biological copies, redundant haplotypes or annotation/assembly redundancy; this assay cannot distinguish those causes. Nonidentical sequences are not proven paralogs. All 526 taxa retained in summary; zero duplicated calls is not evidence against other genome duplication.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'sources'}, indent=2))


if __name__ == '__main__':
    main()
