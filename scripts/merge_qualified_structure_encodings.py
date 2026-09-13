#!/usr/bin/env python3
"""Join unchanged qualified encodings after verifying a disjoint remapped union."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT, sha


def rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def fingerprints(path):
    return Counter(hashlib.sha256(json.dumps(r, sort_keys=True).encode()).digest() for r in rows(path))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', action='append', required=True, help='Mapping directory=qualified encoding directory')
    p.add_argument('--snapshot', type=Path, required=True)
    p.add_argument('--inventory', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable output directory')
    merged = checked_receipt(a.snapshot)
    inventory = checked_receipt(a.inventory)
    if inventory['status'] != 'complete_disjoint_same_method_inventory_union':
        raise ValueError('Expected verified same-method inventory union')
    if merged['source_inventory_sha256'] != sha(a.inventory / 'inventory.jsonl'):
        raise ValueError('Combined mapping uses a different inventory')
    expected_sources = {s['mapping']: s['mapping_receipt_sha256'] for s in inventory['source_cohorts']}
    seen_sources = set(); models = {}; sequences = set(); sources = []; totals = Counter()
    union_tables = {name: Counter() for name in ['marker_structure_links.tsv', 'matrix_to_structure_residues.tsv.gz']}
    provenance = {}
    for arg in a.cohort:
        mapping_name, encoding_name = arg.split('=', 1)
        mapping, encoding = Path(mapping_name), Path(encoding_name)
        mr, er = checked_receipt(mapping), checked_receipt(encoding)
        digest = sha(mapping / 'receipt.json')
        if mapping_name in seen_sources or expected_sources.get(mapping_name) != digest:
            raise ValueError('Unexpected or duplicate source mapping')
        seen_sources.add(mapping_name)
        if er['mapping_receipt_sha256'] != digest or er['confidence_stage'] != 'plddt_and_pae':
            raise ValueError('Qualification does not match mapping or lacks PAE')
        if mr['matrix_receipt_sha256'] != merged['matrix_receipt_sha256']:
            raise ValueError('Different source alignment matrix')
        source_prov = json.loads((mapping / 'model_provenance.json').read_text())
        source_models = {Path(m['path']).stem: m for m in source_prov}
        found = set()
        for row in rows(encoding / 'model_summary.tsv'):
            name = row['model_name']; sid = row['sequence_sha256']
            if name in models or sid in sequences or name not in source_models or name in found:
                raise ValueError('Overlapping or unexpected model/sequence')
            if sid != source_models[name]['sequence_sha256'] or sha(ROOT / row['encoding_path']) != row['encoding_sha256']:
                raise ValueError('Changed encoding or sequence identity')
            models[name] = row; sequences.add(sid); found.add(name)
            for key in er['totals']:
                totals[key] += int(row[key])
        if found != set(source_models) or len(found) != er['models']:
            raise ValueError('Source model grid differs')
        for m in source_prov:
            if m['model_id'] in provenance:
                raise ValueError('Duplicate provenance')
            provenance[m['model_id']] = m
        for filename in union_tables:
            union_tables[filename].update(fingerprints(mapping / filename))
        sources.append({'mapping': mapping_name, 'mapping_receipt_sha256': digest,
                        'encodings': encoding_name, 'encoding_receipt_sha256': sha(encoding / 'receipt.json'),
                        'models': len(found), 'native_config_sha256': er['native_config_sha256'],
                        'coordinate_audit_receipt_sha256': er['coordinate_audit_receipt_sha256'],
                        'pae_receipt_sha256': er['pae_receipt_sha256']})
    if seen_sources != set(expected_sources):
        raise ValueError('Missing source cohort')
    combined_prov = json.loads((a.snapshot / 'model_provenance.json').read_text())
    if len(combined_prov) != len(provenance) or {m['model_id']: m for m in combined_prov} != provenance:
        raise ValueError('Combined model provenance differs from source union')
    counts = {}
    for filename, expected in union_tables.items():
        actual = fingerprints(a.snapshot / filename)
        if actual != expected or any(n != 1 for n in actual.values()):
            raise ValueError('Combined mapping table differs from disjoint source union: ' + filename)
        counts[filename] = sum(actual.values())
    if len(models) != merged['distinct_models'] or len(models) != inventory['models']:
        raise ValueError('Combined model count differs')
    a.output.mkdir(parents=True)
    target = a.output / 'model_summary.tsv'
    with target.open('w') as handle:
        writer = csv.DictWriter(handle, list(next(iter(models.values()))), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(models[k] for k in sorted(models))
    receipt = {'status': 'complete_union_of_previously_qualified_encodings', 'models': len(models),
               'totals': dict(totals), 'confidence_stage': 'plddt_and_pae',
               'mapping_receipt_sha256': sha(a.snapshot / 'receipt.json'),
               'inventory_receipt_sha256': sha(a.inventory / 'receipt.json'),
               'source_cohorts': sources, 'verified_mapping_union_rows': counts,
               'script_sha256': sha(Path(__file__)), 'artifacts': {'model_summary.tsv': sha(target)},
               'interpretation': 'Unchanged encodings reference original audited arrays and retain per-cohort native, coordinate, PAE and qualification receipts. Every combined mapping/provenance row equals the disjoint source union. This is a derived integration receipt, not a new coordinate or PAE audit. Shared settings do not exclude prediction batch effects.'}
    (a.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == '__main__':
    main()
