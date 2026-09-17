#!/usr/bin/env python3
"""Check every completed new-family tree with OrthoFinder's native tree parser."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from orthofinder.tools import tree as native_tree


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    plan = json.loads(a.plan.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    root = Path(plan['output'])
    state = json.loads((root / 'state.json').read_text())
    assert state['status'] == 'complete_family_jobs_requires_independent_readback'
    assert not state['failures'] and state['completed'] == plan['expected_families']
    with (Path(plan['inputs']) / 'family_sequences.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == plan['expected_families']
    assert {p.name for p in (root / 'families').iterdir()} == {r['membership_sha256'] for r in rows}
    tips_total = edges_total = 0
    receipts = {}
    for row in rows:
        key = row['membership_sha256']
        folder = root / 'families' / key
        receipt = json.loads((folder / 'receipt.json').read_text())
        assert receipt['plan_sha256'] == sha(a.plan)
        assert receipt['source_sha256'] == row['sha256']
        source = Path(plan['inputs']) / row['relative_path']
        assert sha(source) == row['sha256']
        for name, digest in receipt['artifacts'].items():
            assert sha(folder / name) == digest
        with source.open() as handle:
            names = [line[1:].strip() for line in handle if line.startswith('>')]
        assert len(set(names)) == len(names) == int(row['proteins'])
        text = (folder / 'tree.nwk').read_text().strip()
        assert text.count(';') == 1 and text.endswith(';')
        tree = native_tree.Tree(text, format=0)
        tips = tree.get_leaf_names()
        assert len(tips) == len(set(tips)) and set(tips) == set(names)
        edges = 0
        for node in tree.traverse():
            if node.is_root():
                continue
            assert math.isfinite(node.dist) and node.dist >= 0
            edges += 1
        assert text.count(':') in [edges, edges + 1]
        tips_total += len(tips)
        edges_total += edges
        receipts[key] = sha(folder / 'receipt.json')
    table = a.output.with_name('native_readback_family_receipts.json')
    if table.exists():
        raise FileExistsError(table)
    table.write_text(json.dumps(receipts, sort_keys=True) + '\n')
    result = dict(status='passed_complete_native_new_family_tree_readback', families=len(rows),
                  tips=tips_total, edges=edges_total, plan_sha256=sha(a.plan),
                  script_sha256=sha(Path(__file__)), artifacts={table.name: sha(table)},
                  scope='Every tree checked with native parser for exact source tips and explicit finite nonnegative lengths; all completed alignment/tree artifacts rehashed. Does not establish biological accuracy, reconciliation or final gene-tree support.')
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
