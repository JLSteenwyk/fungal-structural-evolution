#!/usr/bin/env python3
"""Select a reproducible longest-protein baseline within uniquely mapped genes."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]


def choose(rows):
    groups = defaultdict(list)
    selected, decisions = set(), {}
    for row in rows:
        pid = row['protein_id']
        if pid in decisions:
            raise ValueError('Duplicate protein mapping')
        decisions[pid] = 'pending'
        genes = json.loads(row['gene_ids_json'])
        if row['status'] == 'unique_gene':
            if len(genes) != 1:
                raise ValueError('Inconsistent unique gene mapping')
            groups[genes[0]].append(row)
        else:
            # Preserve unresolved proteins; never infer a shared gene from sequence identity.
            selected.add(pid)
            decisions[pid] = 'retained_unresolved_gene'
    for gene_rows in groups.values():
        ranked = sorted(gene_rows, key=lambda r: (-int(r['protein_length']), r['protein_id']))
        winner = ranked[0]['protein_id']
        selected.add(winner)
        decisions[winner] = 'longest_per_gene_lexical_tiebreak'
        for row in ranked[1:]:
            decisions[row['protein_id']] = 'alternative_product_retained_in_source'
    return selected, decisions


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new output directory for each input snapshot')
    snapshot = json.loads((ROOT / 'metadata/gene_mapping_snapshot.json').read_text())
    inputs = {r['taxon_id']: r for r in json.loads((ROOT / 'metadata/qc_input_receipts.json').read_text())}
    args.output.mkdir(parents=True)
    results = []
    for taxon in snapshot['taxa']:
        name = taxon['taxon_id']
        mapping_path = ROOT / taxon['mapping_path']
        source = ROOT / inputs[name]['input_path']
        if sha(mapping_path) != taxon['mapping_sha256'] or sha(source) != taxon['proteome_sha256']:
            raise ValueError(f'Changed mapping or proteome: {name}')
        with mapping_path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        selected, decisions = choose(rows)
        with source.open() as handle:
            records = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        if set(records) != set(decisions):
            raise ValueError(f'Mapping/proteome accession mismatch: {name}')
        if any(len(records[r['protein_id']].seq) != int(r['protein_length']) for r in rows):
            raise ValueError(f'Mapping/proteome length mismatch: {name}')
        target = args.output / f'{name}.faa'
        with target.open('w') as handle:
            SeqIO.write((records[pid] for pid in sorted(selected)), handle, 'fasta')
        decision_path = args.output / f'{name}.decisions.tsv'
        with decision_path.open('w') as handle:
            writer = csv.DictWriter(handle, list(rows[0]) + ['decision'], delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(dict(row, decision=decisions[row['protein_id']]) for row in rows)
        results.append({'taxon_id': name, 'source_proteins': len(records), 'selected_proteins': len(selected),
                        'unresolved_retained': sum(v == 'retained_unresolved_gene' for v in decisions.values()),
                        'source_sha256': taxon['proteome_sha256'], 'mapping_sha256': taxon['mapping_sha256'],
                        'path': str(target), 'sha256': sha(target), 'decisions_sha256': sha(decision_path)})
    receipt = {'policy': 'longest_protein_per_unique_gene_lexical_accession_ties_v1',
               'status': 'complete' if len(results) == snapshot['planned_taxa'] else 'incomplete_staging',
               'planned_taxa': snapshot['planned_taxa'], 'processed_taxa': len(results), 'taxa': results,
               'mapping_snapshot_sha256': sha(ROOT / 'metadata/gene_mapping_snapshot.json'),
               'note': 'Representative is a baseline, not a canonical isoform claim. Unresolved gene mappings cannot support duplication counts.'}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print('Taxa:', len(results), 'selected:', sum(r['selected_proteins'] for r in results),
          'alternative products:', sum(r['source_proteins'] - r['selected_proteins'] for r in results))


if __name__ == '__main__':
    main()
