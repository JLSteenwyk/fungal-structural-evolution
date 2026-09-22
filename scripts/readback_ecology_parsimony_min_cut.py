#!/usr/bin/env python3
"""Check ML and leave-one-out ecological scores using graph minimum cuts."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import networkx as nx
from Bio import Phylo


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    receipt_path = Path(plan['output'])/'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['plan_sha256'] == sha(args.plan)
    for path, digest in plan['pins'].items():
        assert sha(path) == digest
    table = Path(plan['output'])/'minimum_changes.tsv'
    assert sha(table) == receipt['artifacts'][table.name]
    with open(plan['evidence']) as handle:
        evidence = list(csv.DictReader(handle, delimiter='\t'))
    states = {r['taxon_id']: int(r['state']=='ectomycorrhizal') for r in evidence
              if r['state'] in ('ectomycorrhizal', 'asymbiotic', 'saprotrophic')}
    with table.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    count = 0
    for source in plan['trees']:
        tree = Phylo.read(source['tree'], 'newick')
        nodes = list(tree.find_clades())
        ids = {node:i for i,node in enumerate(nodes)}
        selected = [r for r in rows if r['tree_source']==source['label'] and r['kind']!='ultrafast_bootstrap']
        assert len(selected) == len(states)+1
        assert {r['omitted_taxon'] for r in selected} == set(states)|{''}
        for row in selected:
            graph = nx.Graph()
            for parent in nodes:
                for child in parent.clades:
                    graph.add_edge(ids[parent], ids[child], capacity=1)
            for node in tree.get_terminals():
                if node.name in states and node.name != row['omitted_taxon']:
                    graph.add_edge('state_'+str(states[node.name]), ids[node], capacity=10000)
            value, _ = nx.minimum_cut(graph, 'state_0', 'state_1')
            assert value == int(row['minimum_changes'])
            count += 1
    result = dict(status='passed_independent_graph_min_cut_readback', checked_rows=count,
                  scope='Both maximum-likelihood trees and all leave-one-classification-unknown cases; bootstrap scores not independently recomputed.',
                  algorithm='Unit-capacity tree edges; 10000-capacity terminal class constraints; minimum source-sink cut.',
                  networkx_version=nx.__version__, source_receipt_sha256=sha(receipt_path),
                  script_sha256=sha(__file__))
    args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
