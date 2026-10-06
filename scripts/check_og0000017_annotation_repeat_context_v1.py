#!/usr/bin/env python3
"""Independently check the complete OG0000017 annotation-context audit output."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


FIELDS = ['protein_id', 'native_id', 'group', 'sequence_sha256', 'protein_length', 'gene_id',
          'contig', 'contig_length', 'start', 'end', 'span', 'coding_bases', 'segments',
          'hypothetical', 'scaffold_gene_count', 'left_gap', 'right_gap', 'overlaps_neighbor',
          'uniprot_match', 'pfam_hit_count', 'repeat_keyword_pfam']


def digest(path):
    current = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            current.update(block)
    return current.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--table', type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads(args.receipt.read_text())
    if receipt['status'] != 'complete_og0000017_source_annotation_repeat_risk_audit':
        raise ValueError('Unexpected receipt status')
    if receipt['protein_context_table_sha256'] != digest(args.table):
        raise ValueError('Context-table hash mismatch')
    if receipt['protein_context_table'] != str(args.table):
        raise ValueError('Receipt table path mismatch')
    seen_proteins, seen_native = set(), set()
    groups = {'family': 0, 'other': 0}
    with args.table.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != FIELDS:
            raise ValueError('Unexpected context-table schema')
        for row in reader:
            if row['group'] not in groups:
                raise ValueError('Unexpected cohort')
            if row['protein_id'] in seen_proteins or row['native_id'] in seen_native:
                raise ValueError('Repeated protein or native identifier')
            if not row['native_id'].startswith('214_'):
                raise ValueError('Non-focal native identifier')
            if len(row['sequence_sha256']) != 64:
                raise ValueError('Invalid sequence hash')
            for numeric in ('protein_length', 'contig_length', 'start', 'end', 'span',
                            'coding_bases', 'segments', 'scaffold_gene_count', 'pfam_hit_count'):
                if int(row[numeric]) < 0:
                    raise ValueError('Negative numeric value')
            seen_proteins.add(row['protein_id']); seen_native.add(row['native_id'])
            groups[row['group']] += 1
    if len(seen_proteins) != receipt['source_products'] or len(seen_proteins) != 133973:
        raise ValueError('Source product census mismatch')
    if groups['family'] != receipt['focal_family_products'] or groups['family'] != 46787:
        raise ValueError('Family census mismatch')
    if groups['other'] != receipt['other_focal_products'] or groups['other'] != 87186:
        raise ValueError('Other-product census mismatch')
    print(json.dumps({'status': 'passed_og0000017_annotation_repeat_context_output_check',
                      'source_products': len(seen_proteins), 'cohorts': groups,
                      'table_sha256': receipt['protein_context_table_sha256']}, indent=2))


if __name__ == '__main__':
    main()
