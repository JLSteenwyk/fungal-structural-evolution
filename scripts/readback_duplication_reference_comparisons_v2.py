#!/usr/bin/env python3
"""Check the complete reference ledger against native genes and frozen model records."""
import argparse
import csv
import hashlib
import json
import sqlite3
from collections import Counter
from pathlib import Path

from run_ortholog_pair_guide_comparison import sha


def identity(left, right):
    endpoints = tuple(sorted([left, right]))
    return endpoints, hashlib.sha256(json.dumps(endpoints, separators=(',', ':')).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['inventory', 'references', 'reference-readback', 'base-queue',
                 'base-readback', 'bridge', 'catalog', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    receipt_path = args.inventory / 'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    if receipt['status'] != 'complete_provisional_reference_comparison_inventory':
        raise ValueError('Incomplete comparison inventory')
    references_path = args.references / 'receipt.json'
    references = json.loads(references_path.read_text())
    reference_proof = json.loads(args.reference_readback.read_text())
    if (references['status'] != 'complete_duplication_sister_reference_inventory'
            or reference_proof['status'] != 'passed_full_duplication_sister_reference_readback'
            or reference_proof['producer_receipt_sha256'] != sha(references_path)):
        raise ValueError('Unverified reference choices')
    base_path = args.base_queue / 'receipt.json'
    base = json.loads(base_path.read_text())
    base_proof = json.loads(args.base_readback.read_text())
    if (base['status'] != 'complete_reviewed_duplication_model_pair_queue'
            or base_proof['status'] != 'passed_full_duplication_model_pair_queue_export_readback'):
        raise ValueError('Unverified primary queue')
    bindings = dict(receipt['source_pins'])
    for path in [references_path, base_path, args.bridge, args.catalog / 'receipt.json']:
        if bindings.get(str(path)) != sha(path):
            raise ValueError('CLI source differs from producer: ' + str(path))
    for path in [receipt_path, args.reference_readback, args.base_readback]:
        bindings[str(path)] = sha(path)
    for root, record in [(args.inventory, receipt), (args.references, references)]:
        for name, digest in record['artifacts'].items():
            path = str(root / name)
            if path in bindings and bindings[path] != digest:
                raise ValueError('Conflicting artifact binding')
            bindings[path] = digest
    for name, digest in reference_proof.get('artifacts', {}).items():
        bindings[str(args.reference_readback.parent / name)] = digest
    for name in ['receipt.json', 'models.jsonl', 'model_pairs.tsv', 'event_model_pair_links.tsv']:
        path = str(args.base_queue / name)
        if base_proof['source_sha256'].get(path) != sha(path):
            raise ValueError('Primary queue/readback binding differs: ' + name)
        bindings[path] = sha(path)
        if name != 'receipt.json' and bindings[path] != base['artifacts'][name]:
            raise ValueError('Changed primary queue artifact')
    catalog = json.loads((args.catalog / 'receipt.json').read_text())
    catalog_path = args.catalog / 'models.jsonl'
    if bindings.get(str(catalog_path)) != catalog['artifacts']['models.jsonl']:
        raise ValueError('Catalog not bound to comparison producer')

    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed pinned source or output: ' + path)

    verify()
    expected = {}
    source_counts = Counter()
    needed = set()
    for guide in ['profile', 'mafft']:
        source_rows = 0
        with (args.references / (guide + '_sister_references.tsv')).open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                source_rows += 1
                if row['status'] != 'provisional_reference_available':
                    continue
                ties = json.loads(row['nearest_reference_genes'])
                if not ties or ties != sorted(set(ties)) or row['chosen_reference_gene'] != ties[0]:
                    raise ValueError('Invalid reference ties')
                source_counts[guide + '_eligible_events'] += 1
                source_counts[guide + '_event_reference_links'] += len(ties)
                for reference in ties:
                    if reference.split('_', 1)[0] == row['taxon_id']:
                        raise ValueError('Focal taxon used as reference')
                    for side in ['a', 'b']:
                        key = (guide, row['family'], row['gene_node'], row['gene_a'], row['gene_b'],
                               reference, side, row['gene_' + side], str(int(reference == ties[0])))
                        if key in expected:
                            raise ValueError('Repeated source event/reference link')
                        expected[key] = (row['gene_' + side], reference)
                        needed.update(expected[key])
        summary = next(item for item in references['guides'] if item['guide'] == guide)
        audited = next(item for item in reference_proof['guides'] if item['guide'] == guide)
        if source_rows != summary['candidates'] or summary['candidates'] != audited['candidates']:
            raise ValueError('Reference source universe differs')
    links = {}
    with sqlite3.connect('file:' + str(args.bridge.resolve()) + '?mode=ro', uri=True) as connection:
        for taxon, protein, sequence, model, version, path in connection.execute(
                'SELECT taxon_id,protein_id,sequence_sha256,model_id,version,model_path FROM structures'):
            gene = taxon + '_' + protein
            if gene in needed:
                if gene in links:
                    raise ValueError('Repeated native gene/model mapping')
                links[gene] = (model, version, sequence, path)
    if set(links) != needed:
        raise ValueError('Missing native gene/model mapping')
    base_models = set()
    with (args.base_queue / 'models.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line)
            key = row['model_id'], row['version']
            if key in base_models:
                raise ValueError('Repeated primary model')
            base_models.add(key)
    if len(base_models) != base_proof['models_checked']:
        raise ValueError('Primary model scope differs')
    base_pairs = set()
    with (args.base_queue / 'model_pairs.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            ends, digest = identity((row['model_a'], int(row['version_a'])),
                                    (row['model_b'], int(row['version_b'])))
            if digest != row['pair_key'] or ends in base_pairs or ends[0] == ends[1]:
                raise ValueError('Invalid primary pair')
            base_pairs.add(ends)
    if len(base_pairs) != base_proof['distinct_pairs_checked']:
        raise ValueError('Primary pair scope differs')
    observed = set()
    models = set()
    pairs = {}
    work_counts = Counter()
    fields = ['guide', 'family', 'gene_node', 'gene_a', 'gene_b', 'reference_gene',
              'focal_side', 'focal_gene', 'lexical_representative']
    with (args.inventory / 'event_reference_comparisons.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = tuple(row[name] for name in fields)
            if key in observed or key not in expected:
                raise ValueError('Unknown/repeated event-reference link')
            observed.add(key)
            focal, reference = expected[key]
            left = row['focal_model'], int(row['focal_version'])
            right = row['reference_model'], int(row['reference_version'])
            if left != links[focal][:2] or right != links[reference][:2]:
                raise ValueError('Native gene/model/version assignment differs')
            models.update([left, right])
            ends, digest = identity(left, right)
            status = ('identical_model' if left == right else
                      'existing_duplicate_pair' if ends in base_pairs else 'additional_pair')
            if row['pair_key'] != digest or row['work_disposition'] != status:
                raise ValueError('Pair identity or work partition differs')
            work_counts[status] += 1
            if left != right:
                pairs[digest] = ends
    if observed != set(expected):
        raise ValueError('Complete tied-reference ledger differs')
    expected_pairs = {}
    for digest, (left, right) in pairs.items():
        expected_pairs[digest] = dict(pair_key=digest, model_a=left[0], version_a=str(left[1]),
                                     model_b=right[0], version_b=str(right[1]),
                                     work_disposition='existing_duplicate_pair' if (left, right) in base_pairs else 'additional_pair')
    written = set()
    with (args.inventory / 'model_pairs.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['pair_key']
            if key in written or key not in expected_pairs or row != expected_pairs[key]:
                raise ValueError('Exported complete pair row differs')
            written.add(key)
    if written != set(pairs):
        raise ValueError('Missing pair output')
    selected = {}
    with catalog_path.open() as handle:
        for line in handle:
            row = json.loads(line)
            key = row['model_id'], row['version']
            if key in models:
                if key in selected:
                    raise ValueError('Repeated source catalog identity')
                selected[key] = row
    if set(selected) != models:
        raise ValueError('Missing catalog model')
    for model, version, sequence, path in links.values():
        row = selected[model, version]
        if (row['sequence_sha256'], row['path']) != (sequence, path):
            raise ValueError('Catalog gene/sequence/path binding differs')
    extras = models - base_models
    for name, desired in [('models.jsonl', models), ('additional_models.jsonl', extras)]:
        found = set()
        with (args.inventory / name).open() as handle:
            for line in handle:
                row = json.loads(line)
                key = row['model_id'], row['version']
                if key not in desired or key in found or row != selected[key]:
                    raise ValueError('Exported catalog model record differs')
                found.add(key)
        if found != desired:
            raise ValueError('Incomplete model partition')
    additional_pairs = sum(ends not in base_pairs for ends in pairs.values())
    counts = {**source_counts, **{status + '_event_links': n for status, n in work_counts.items()}}
    byte_count = sum(Path(selected[key]['path']).stat().st_size for key in extras)
    expected_counts = dict(all_reference_comparison_models=len(models), additional_models=len(extras),
                           additional_coordinate_bytes=byte_count, unique_distinct_model_pairs=len(pairs),
                           additional_model_pairs=additional_pairs,
                           existing_duplicate_model_pairs=len(pairs) - additional_pairs,
                           additional_directed_mask_dispositions=4 * additional_pairs)
    if receipt['counts'] != counts or any(receipt[name] != n for name, n in expected_counts.items()):
        raise ValueError('Full source/model/work summary differs')
    verify()
    result = dict(status='passed_full_reference_comparison_ledger_and_native_model_readback',
                  producer_receipt_sha256=sha(receipt_path),
                  reference_readback_sha256=sha(args.reference_readback),
                  primary_queue_readback_sha256=sha(args.base_readback),
                  event_reference_comparisons=len(observed), unique_models=len(models),
                  additional_models=len(extras), distinct_model_pairs=len(pairs),
                  additional_model_pairs=additional_pairs, work_dispositions=dict(work_counts),
                  native_gene_model_mappings_checked=len(links), counts=counts,
                  script_sha256=sha(__file__), source_hashes=bindings,
                  scope='Every event/tied-reference/focal-side identity, native gene/model/version assignment, '
                  'catalog sequence/path and complete exported model record, distinct pair row and work '
                  'partition independently reconstructed. Reference choices require separate full native-tree '
                  'readback. Coordinate byte counts are stat checks, not raw coordinate validation, independent '
                  'biological orthology, ancestral states or structural asymmetry inference.')
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
