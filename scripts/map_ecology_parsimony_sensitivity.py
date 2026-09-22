#!/usr/bin/env python3
"""Unrooted minimum-change diagnostic for published ecological classifications."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from Bio import Phylo


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def minimum_changes(tree, states):
    """Binary symmetric Sankoff cost; unspecified tips allow either state."""
    costs = {}
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal():
            state = states.get(node.name)
            costs[node] = (0, 0) if state is None else ((0, 10**9) if state == 0 else (10**9, 0))
        else:
            costs[node] = tuple(sum(min(costs[c][s], costs[c][1-s]+1)
                                    for c in node.clades) for s in (0, 1))
    return min(costs[tree.root])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    for path, digest in plan['pins'].items():
        if sha(path) != digest:
            raise ValueError('Changed input: '+path)
    out = Path(plan['output'])
    if out.exists():
        raise FileExistsError(out)
    out.mkdir(parents=True)
    with open(plan['evidence']) as handle:
        evidence = list(csv.DictReader(handle, delimiter='\t'))
    with open(plan['manifest']) as handle:
        manifest = list(csv.DictReader(handle, delimiter='\t'))
    taxa = {r['taxon_id'] for r in manifest}
    assert len(taxa) == len(manifest)
    states = {}
    for row in evidence:
        if row['state'] == 'ectomycorrhizal':
            states[row['taxon_id']] = 1
        elif row['state'] in ('asymbiotic', 'saprotrophic'):
            states[row['taxon_id']] = 0
        elif row['state'] != 'orchid_mycorrhizal':
            raise ValueError('Unreviewed state coding')
    assert set(states) <= taxa
    rows = []
    summaries = []
    for source in plan['trees']:
        receipt = json.loads(Path(source['receipt']).read_text())
        audit = json.loads(Path(source['readback']).read_text())
        assert audit['status'] == 'passed_pmsf_profile_tree_and_bootstrap_readback'
        assert audit['source_receipt_sha256'] == sha(source['receipt'])
        for path in (source['tree'], source['bootstrap']):
            assert sha(path) == receipt['artifacts'][Path(path).name]
        def check(tree):
            tips = [n.name for n in tree.get_terminals()]
            assert len(tips) == len(taxa) and set(tips) == taxa
        tree = Phylo.read(source['tree'], 'newick')
        check(tree)
        baseline = minimum_changes(tree, states)
        rows.append(dict(tree_source=source['label'], kind='maximum_likelihood', index=0,
                         omitted_taxon='', minimum_changes=baseline))
        for taxon in sorted(states):
            value = minimum_changes(tree, {t:s for t,s in states.items() if t != taxon})
            assert value <= baseline
            rows.append(dict(tree_source=source['label'], kind='leave_one_classification_unknown',
                             index=0, omitted_taxon=taxon, minimum_changes=value))
        distribution = Counter()
        for number, bootstrap in enumerate(Phylo.parse(source['bootstrap'], 'newick'), 1):
            check(bootstrap)
            value = minimum_changes(bootstrap, states)
            distribution[value] += 1
            rows.append(dict(tree_source=source['label'], kind='ultrafast_bootstrap', index=number,
                             omitted_taxon='', minimum_changes=value))
        assert sum(distribution.values()) == audit['bootstrap_trees'] == 1000
        summaries.append(dict(source=source['label'], ml_minimum_changes=baseline,
                              bootstrap_distribution=dict(sorted(distribution.items()))))
        print(source['label'], summaries[-1], flush=True)
    output = out/'minimum_changes.tsv'
    with output.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    for path, digest in plan['pins'].items():
        assert sha(path) == digest
    result = dict(status='complete_ecological_classification_parsimony_diagnostic',
                  plan_sha256=sha(args.plan), states=dict(Counter(states.values())),
                  unknown_tips=len(taxa)-len(states), summaries=summaries, rows=len(rows),
                  artifacts={'minimum_changes.tsv':sha(output)},
                  scope='Minimum total changes under binary symmetric equal costs, conditional on published classifications. Unknowns remain unconstrained, including orchid symbiosis. No rooted gains/losses, replicated-origin count, transition locations, effect test or statistical confidence interval. Bootstrap distribution samples topology variation on one sequence alignment, not trait or model uncertainty.')
    (out/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
