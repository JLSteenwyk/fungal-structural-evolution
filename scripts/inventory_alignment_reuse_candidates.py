#!/usr/bin/env python3
"""Inventory the full old/new pair union; matching catalog bytes do not authorize reuse."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

from screen_duplication_alignment_reuse import sha


FIELDS = ['pair_key', 'model_a', 'version_a', 'model_b', 'version_b',
          'new_sources', 'matching_old_sources', 'changed_old_sources', 'disposition']


def load_sources(specs):
    pairs, counts = {}, {}
    for spec in specs:
        label = spec['label']
        if label in counts:
            raise ValueError('Repeated source label')
        root = Path(spec['root'])
        receipt = json.loads((root / 'receipt.json').read_text())
        if receipt['status'] != spec['status']:
            raise ValueError('Incomplete source: ' + label)
        models = {}
        for name in ['models.jsonl', 'model_pairs.tsv']:
            if sha(root / name) != receipt['artifacts'][name]:
                raise ValueError('Changed source artifact: ' + label + ':' + name)
        with (root / 'models.jsonl').open() as handle:
            for line in handle:
                row = json.loads(line)
                key = row['model_id'], int(row['version'])
                if key in models or row['version'] != key[1] or row['length'] <= 0:
                    raise ValueError('Invalid model identity')
                models[key] = (row['sha256'], row['sequence_sha256'], row['length'])
        seen = set()
        with (root / 'model_pairs.tsv').open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                ends = tuple(sorted([(row['model_a'], int(row['version_a'])),
                                     (row['model_b'], int(row['version_b']))]))
                digest = hashlib.sha256(json.dumps(ends, separators=(',', ':')).encode()).hexdigest()
                if (digest != row['pair_key'] or digest in seen or ends[0] == ends[1]
                        or any(key not in models for key in ends)):
                    raise ValueError('Invalid pair identity')
                seen.add(digest)
                signature = tuple(models[key] for key in ends)
                record = pairs.setdefault(digest, dict(ends=ends, sources={}))
                if record['ends'] != ends:
                    raise ValueError('Pair hash collision')
                record['sources'][label] = signature
        if len(models) != receipt[spec['model_count_field']] or len(seen) != receipt[spec['pair_count_field']]:
            raise ValueError('Source counts differ')
        counts[label] = dict(models=len(models), pairs=len(seen))
    return pairs, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}

    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed pin: ' + path)

    for spec in plan['old'] + plan['new']:
        for name in ['receipt.json', 'models.jsonl', 'model_pairs.tsv']:
            if str(Path(spec['root']) / name) not in bindings:
                raise ValueError('Unpinned source')
    verify()
    old, old_counts = load_sources(plan['old'])
    new, new_counts = load_sources(plan['new'])
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    counts, per_source = Counter(), Counter()
    with (out / 'pair_reuse_candidates.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for digest, record in sorted(new.items()):
            signatures = list(record['sources'].values())
            if any(signature != signatures[0] for signature in signatures):
                raise ValueError('Conflicting frozen new source descriptors')
            previous = old.get(digest, dict(ends=record['ends'], sources={}))
            if previous['ends'] != record['ends']:
                raise ValueError('Old/new pair endpoints differ')
            matching = sorted(k for k, v in previous['sources'].items() if v == signatures[0])
            changed = sorted(set(previous['sources']) - set(matching))
            status = ('matching_catalog_sources_pending_input_and_result_checks' if matching else
                      'changed_catalog_sources_require_alignment' if changed else
                      'new_pair_requires_alignment')
            (left, lv), (right, rv) = record['ends']
            writer.writerow(dict(pair_key=digest, model_a=left, version_a=lv,
                                 model_b=right, version_b=rv,
                                 new_sources=json.dumps(sorted(record['sources']), separators=(',', ':')),
                                 matching_old_sources=json.dumps(matching, separators=(',', ':')),
                                 changed_old_sources=json.dumps(changed, separators=(',', ':')),
                                 disposition=status))
            counts[status] += 1
            for label in record['sources']:
                per_source[label + ':' + status] += 1
    verify()
    receipt = dict(status='complete_full_pair_union_catalog_reuse_screen_not_authorization',
                   old_union_pairs=len(old), new_union_pairs=len(new),
                   old_union_pairs_absent=len(set(old) - set(new)), counts=dict(counts),
                   per_new_source_counts=dict(per_source), old_sources=old_counts, new_sources=new_counts,
                   source_hashes=bindings, script_sha256=sha(__file__),
                   artifacts={'pair_reuse_candidates.tsv': sha(out / 'pair_reuse_candidates.tsv')},
                   scope='Complete model/version pair union across all declared sources, deduplicating overlapping inventories. Exact catalog coordinate checksum, sequence checksum and length compared per endpoint; every matching/changed previous source retained. No raw coordinate rehash, input/mask/order/executable/checkpoint verification, numeric-flag clearance, result copying or reuse authorization. Native structural inference and quality qualification remain separate stages.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ['status', 'old_union_pairs', 'new_union_pairs', 'counts']}), flush=True)


if __name__ == '__main__':
    main()
