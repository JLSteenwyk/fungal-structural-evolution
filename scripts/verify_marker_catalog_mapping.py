#!/usr/bin/env python3
"""Verify exact agreement of the prefetch catalog with a completed residue mapping."""
import argparse
import csv
import json
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha


def keyed(rows, keys, fields):
    result = {}
    for row in rows:
        key = tuple(row[k] for k in keys)
        if key in result:
            raise ValueError('Duplicate identity')
        result[key] = tuple(str(row[k]) for k in fields)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog', type=Path, required=True)
    parser.add_argument('--mapping', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new verification receipt')
    catalog = checked_receipt(args.catalog)
    mapping = checked_receipt(args.mapping)
    if catalog['inventory_sha256'] != mapping['source_inventory_sha256']:
        raise ValueError('Different inventory snapshots')
    if any(catalog['source_policy'][k] != mapping['source_policy'][k] for k in ['provider', 'tool']):
        raise ValueError('Different source policies')
    models = [json.loads((p / 'model_provenance.json').read_text()) for p in [args.catalog, args.mapping]]
    fields = ['path', 'sha256', 'sequence_sha256', 'length', 'provider', 'tool']
    model_keys = [keyed(rows, ['model_id', 'version'], fields) for rows in models]
    if model_keys[0] != model_keys[1]:
        raise ValueError('Model identities or coordinate provenance differ')
    with (args.catalog / 'marker_model_links.tsv').open() as handle:
        original = list(csv.DictReader(handle, delimiter='\t'))
    with (args.mapping / 'marker_structure_links.tsv').open() as handle:
        final = [dict(r, version=r['model_version']) for r in csv.DictReader(handle, delimiter='\t')]
    fields = ['protein_id', 'sequence_sha256', 'model_id', 'version', 'model_path']
    if keyed(original, ['marker', 'taxon_id'], fields) != keyed(final, ['marker', 'taxon_id'], fields):
        raise ValueError('Marker linkage differs')
    if len(models[0]) != mapping['distinct_models'] or len(final) != mapping['marker_proteins_linked']:
        raise ValueError('Receipt count mismatch')
    result = {'status': 'complete_catalog_mapping_identity_agreement', 'models': len(models[0]),
        'marker_links': len(final), 'catalog_path': str(args.catalog), 'mapping_path': str(args.mapping),
        'catalog_receipt_sha256': sha(args.catalog / 'receipt.json'),
        'mapping_receipt_sha256': sha(args.mapping / 'receipt.json'),
        'script_sha256': sha(Path(__file__)),
        'interpretation': 'All receipt artifacts checked; exact model/version, coordinate path/hash, sequence identity, source policy and marker-link agreement. This does not validate PAE or native features. Prefetch PAE must still be cache-validated into a receipt tied to the final mapping before downstream confidence use.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
