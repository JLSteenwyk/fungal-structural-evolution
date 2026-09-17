#!/usr/bin/env python3
"""Validate retained trees and map exact memberships into both expanded guides.

Read-only with respect to source trees. Missing trees remain explicit pending
rows; this catalog does not perform reconciliation or approve biological trees.
Run with the project's native OrthoFinder Python environment.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from orthofinder.tools import tree as native_tree


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['census', 'staged', 'merged', 'output']:
        ap.add_argument('--' + name, type=Path, required=True)
    a = ap.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    census_receipt = json.loads((a.census / 'receipt.json').read_text())
    for name, digest in census_receipt['artifacts'].items():
        if sha(a.census / name) != digest:
            raise ValueError('Changed census artifact: ' + name)
    merge_receipt = json.loads((a.merged / 'receipt.json').read_text())
    if sha(a.merged / 'receipt.json') != census_receipt['merge_receipt_sha256']:
        raise ValueError('Changed merged partition receipt')
    tasks = {r['membership_sha256']: r for r in rows(a.census / 'unique_tree_tasks.tsv')
             if r['source_type'] == 'retained'}
    copied = {r['relative_path']: r for r in rows(a.staged / 'copied_files.tsv')}
    mappings = {key: {} for key in tasks}
    for guide in merge_receipt['guides']:
        for name, digest in guide['artifacts'].items():
            if sha(a.merged / name) != digest:
                raise ValueError('Changed merged artifact: ' + name)
        for row in rows(a.merged / guide['guide'] / 'family_sources.tsv'):
            if row['source_type'] != 'retained' or int(row['proteins']) < 3:
                continue
            key = row['membership_sha256']
            task = tasks[key]
            if (row['original_tree_candidate'] != task['original_tree_candidate']
                    or int(row['proteins']) != int(task['proteins'])
                    or guide['guide'] in mappings[key]):
                raise ValueError('Retained membership crosswalk mismatch')
            mappings[key][guide['guide']] = row['new_family']
    fields = ['membership_sha256', 'original_family', 'profile_family',
              'mafft_family', 'proteins', 'taxa', 'tree_path', 'tree_sha256', 'status']
    passed = 0
    pending = []
    table = a.output / 'retained_tree_catalog.tsv'
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for key, task in tasks.items():
            if set(mappings[key]) != {'profile', 'mafft'}:
                raise ValueError('Missing guide mapping')
            family = task['original_tree_candidate']
            relative = f'Source/WorkingDirectory/Trees_ids/{family}.txt'
            path = a.staged / relative
            record = dict(membership_sha256=key, original_family=family,
                          profile_family=mappings[key]['profile'],
                          mafft_family=mappings[key]['mafft'], proteins=task['proteins'],
                          taxa=task['taxa'], tree_path=str(path.resolve()), tree_sha256='')
            if relative not in copied:
                if path.exists():
                    raise ValueError('Unmanifested staged tree: ' + family)
                record['status'] = 'pending_missing_staged_tree'
                pending.append(family)
            else:
                digest = sha(path)
                if digest != copied[relative]['sha256'] or path.stat().st_size != int(copied[relative]['bytes']):
                    raise ValueError('Changed staged tree: ' + family)
                text = path.read_text().strip()
                if text.count(';') != 1 or not text.endswith(';'):
                    raise ValueError('Incomplete tree: ' + family)
                parsed = native_tree.Tree(text, format=0)
                tips = parsed.get_leaf_names()
                membership = hashlib.sha256(('\n'.join(sorted(tips)) + '\n').encode()).hexdigest()
                if (len(tips) != len(set(tips)) or membership != key
                        or len(tips) != int(task['proteins'])
                        or len({tip.split('_', 1)[0] for tip in tips}) != int(task['taxa'])):
                    raise ValueError('Native tip membership mismatch: ' + family)
                edges = 0
                for node in parsed.traverse():
                    if node.is_root():
                        continue
                    edges += 1
                    if not math.isfinite(node.dist) or node.dist < 0:
                        raise ValueError('Invalid branch length: ' + family)
                if any(':' in tip for tip in tips) or text.count(':') not in [edges, edges + 1]:
                    raise ValueError('Missing explicit edge lengths: ' + family)
                if sha(path) != digest:
                    raise ValueError('Tree changed during validation')
                record.update(tree_sha256=digest, status='validated_exact_membership')
                passed += 1
            writer.writerow(record)
            if (passed + len(pending)) % 2000 == 0:
                print(json.dumps(dict(validated=passed, pending=len(pending))), flush=True)
    result = dict(status='completed_retained_tree_catalog', validated_trees=passed,
                  pending_families=pending, retained_tasks=len(tasks),
                  guide_mappings=2 * len(tasks), catalog_sha256=sha(table),
                  census_receipt_sha256=sha(a.census / 'receipt.json'),
                  merge_receipt_sha256=sha(a.merged / 'receipt.json'),
                  staged_manifest_sha256=sha(a.staged / 'copied_files.tsv'),
                  script_sha256=sha(Path(__file__)), parser_sha256=sha(Path(native_tree.__file__)),
                  scope='Every retained tree candidate mapped to both expanded partitions; available immutable staged trees verified using native parser, exact membership fingerprint and explicit finite nonnegative edges. Pending trees are not validated. No reconciliation or biological topology validation.')
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
