#!/usr/bin/env python3
"""Check repaired FastTree output using OrthoFinder's bundled ETE parser."""
import argparse
import hashlib
import json
import math
from pathlib import Path
from orthofinder.tools import tree as native_tree
from orthofinder.tools import newick


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['tree', 'source', 'output']:
        ap.add_argument('--' + name, type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    source_hash, tree_hash = sha(a.source), sha(a.tree)
    with a.source.open() as handle:
        names = [line[1:].strip() for line in handle if line.startswith('>')]
    if not names or len(names) != len(set(names)):
        raise ValueError('Empty or duplicated source tip universe')
    text = a.tree.read_text().strip()
    if text.count(';') != 1 or not text.endswith(';'):
        raise ValueError('Exactly one complete Newick tree required')
    parsed = native_tree.Tree(text, format=0)
    tips = parsed.get_leaf_names()
    if len(tips) != len(set(tips)) or set(tips) != set(names):
        raise ValueError('Native parser tip universe mismatch')
    lengths = []
    for node in parsed.traverse():
        if node.is_root():
            continue
        if not math.isfinite(node.dist) or node.dist < 0:
            raise ValueError('Nonfinite or negative branch length')
        lengths.append(node.dist)
    # The native reader defaults missing lengths; require every non-root edge
    # to have an explicit serialized length. Source identifiers contain no colons.
    if any(':' in name for name in names) or text.count(':') not in [len(lengths), len(lengths) + 1]:
        raise ValueError('Branch-length serialization incomplete or ambiguous')
    if source_hash != sha(a.source) or tree_hash != sha(a.tree):
        raise ValueError('Source changed during validation')
    result = dict(status='passed_native_repaired_tree_readback', tips=len(tips),
                  edges=len(lengths), total_branch_length=math.fsum(lengths),
                  source_sha256=source_hash, tree_sha256=tree_hash,
                  parser_sha256={str(Path(p)): sha(Path(p)) for p in [native_tree.__file__, newick.__file__]},
                  script_sha256=sha(Path(__file__)),
                  scope='Independent native-parser tip and explicit finite nonnegative branch checks. Does not validate biological topology, orthology or reconciliation.')
    a.output.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
