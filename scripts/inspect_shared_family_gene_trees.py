#!/usr/bin/env python3
"""Inspect focal candidate membership and local context in audited resolved trees."""
import argparse
import csv
from io import StringIO
import json
from pathlib import Path
from Bio import Phylo
from map_cross_clan_families import sha


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pairs', type=Path, required=True)
    ap.add_argument('--profile', type=Path, required=True)
    ap.add_argument('--mafft', type=Path, required=True)
    ap.add_argument('--profile-audit', type=Path, required=True)
    ap.add_argument('--mafft-audit', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    paths = [a.pairs, a.profile, a.mafft, a.profile_audit, a.mafft_audit,
             a.pairs.parent / 'receipt.json', Path(__file__), Path(__file__).with_name('map_cross_clan_families.py')]
    pins = {str(p): sha(p) for p in paths}
    assert json.loads(paths[5].read_text())['artifacts'][a.pairs.name] == pins[str(a.pairs)]
    with a.pairs.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    results = []
    for guide in ['profile', 'mafft']:
        source = getattr(a, guide); audit = json.loads(getattr(a, guide + '_audit').read_text())
        assert audit['status'] == 'passed_complete_resolved_tree_membership_readback'
        assert audit['resolved_tree_file_sha256'] == pins[str(source)]
        selected = [r for r in rows if r['guide'] == guide]
        needed = {r['family'] for r in selected}; trees = {}
        with source.open() as f:
            for line in f:
                family, newick = line.rstrip().split(': ', 1)
                if family in needed:
                    assert family not in trees
                    trees[family] = newick
        assert set(trees) == needed
        for family, newick in sorted(trees.items()):
            tree = Phylo.read(StringIO(newick), 'newick')
            tips = tree.get_terminals(); labels = {t.name: t for t in tips}
            assert len(labels) == len(tips)
            members = [r for r in selected if r['family'] == family]
            sides = {side: {r['taxon_id_' + side] + '_' + r['protein_id_' + side] for r in members} for side in ['left', 'right']}
            assert (sides['left'] | sides['right']) <= set(labels)
            parents = {child: node for node in tree.find_clades() for child in node.clades}
            focal = []
            for label in sorted(sides['right']):
                leaf = labels[label]; parent = parents[leaf]
                sister = sorted(t.name for child in parent.clades if child is not leaf for t in child.get_terminals())
                focal.append(dict(label=label, terminal_branch_length=leaf.branch_length,
                                  parent_label=parent.name, sister_labels=sister,
                                  sister_candidate_left_labels=sorted(set(sister) & sides['left'])))
            mrca = tree.common_ancestor([labels[l] for l in sorted(sides['right'])])
            results.append(dict(guide=guide, family=family, tips=len(tips),
                                candidate_left_tips=len(sides['left']), candidate_right_tips=len(sides['right']),
                                right_mrca_tips=len(mrca.get_terminals()), focal_context=focal))
    for p, digest in pins.items():
        assert sha(p) == digest
    out = dict(status='complete_focal_resolved_tree_membership_context', source_hashes=pins, results=results,
               scope='Exact focal candidate membership and immediate sister context in existing native resolved trees. Branch lengths are sequence divergence, not time or structural displacement. Native node names are not support values. Does not validate rooting, reconciliation events, orthology, domain histories, or structural acceleration.')
    with a.output.open('x') as f:
        json.dump(out, f, indent=2); f.write('\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
