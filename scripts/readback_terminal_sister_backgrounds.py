#!/usr/bin/env python3
"""Reconstruct every terminal sister record independently from native trees."""
import argparse
import csv
import json
import math
import sqlite3
import time
from collections import Counter, defaultdict
from io import StringIO
from pathlib import Path
import psutil
from Bio import Phylo
from run_ortholog_pair_guide_comparison import sha


def reconstruct(newick, reports, lookup, taxa):
    tree = Phylo.read(StringIO(newick), 'newick')
    parents = {child: parent for parent in tree.find_clades() for child in parent.clades}
    seen = set(); result = {}
    for node in tree.find_clades(order='postorder'):
        if not node.name or node.name in seen:
            raise ValueError('Missing or repeated node name')
        seen.add(node.name)
        if node is not tree.root and (node.branch_length is None or not math.isfinite(node.branch_length) or node.branch_length < 0):
            raise ValueError('Invalid branch')
        if len(node.clades) != 2 or not all(child.is_terminal() for child in node.clades):
            continue
        leaves = sorted(node.clades, key=lambda child: child.name)
        genes = [child.name for child in leaves]
        species = [gene.partition('_')[0] for gene in genes]
        assert set(species) <= taxa
        parent = parents.get(node)
        reported = node.name in reports
        label = 'cross_taxon_unreported_candidate'
        if species[0] == species[1]:
            label = 'same_taxon_unreported'
        if reported:
            label = 'reported_duplication'
        models = [lookup.get(gene) for gene in genes]
        n_models = sum(model is not None for model in models)
        coverage = ['no_models', 'one_model', 'two_distinct_models'][n_models]
        if n_models == 2 and models[0] == models[1]:
            coverage = 'identical_model'
        row = dict(gene_node=node.name, parent_node=parent.name if parent else '',
                   node_reported_duplication=str(int(reported)),
                   parent_reported_duplication=str(int(parent is not None and parent.name in reports)),
                   gene_a=genes[0], gene_b=genes[1], taxon_a=species[0], taxon_b=species[1],
                   sequence_tip_a=leaves[0].branch_length, sequence_tip_b=leaves[1].branch_length,
                   sequence_pair_distance=sum(child.branch_length for child in leaves),
                   candidate_status=label, model_coverage=coverage)
        for suffix, model in zip(['a', 'b'], models):
            row['model_' + suffix] = str(model[0]) if model else ''
            row['version_' + suffix] = str(model[1]) if model else ''
        result[node.name] = row
    return result


def compare(actual, expected):
    assert actual.keys() == expected.keys(), 'Field set differs'
    for field, value in expected.items():
        if field.startswith('sequence_'):
            assert float(actual[field]) == value, field
        else:
            assert actual[field] == value, field


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text()); config_hash = sha(args.config)
    def verify_config():
        assert sha(args.config) == config_hash
        for path, digest in config['pins'].items():
            assert sha(path) == digest, path
    verify_config()
    dep = config['dependency']
    while True:
        try:
            process = psutil.Process(dep['pid'])
            if process.create_time() != dep['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            assert process.cmdline() == dep['cmdline'], 'Producer command changed'
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    verify_config()
    plan = json.loads(Path(config['plan']).read_text())
    root = Path(plan['output']); receipt_path = root / 'receipt.json'
    receipt_hash = sha(receipt_path); receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'complete_terminal_sister_inventory_pending_independent_readback'
    assert receipt['plan_sha256'] == sha(config['plan'])
    def verify_inputs():
        assert sha(receipt_path) == receipt_hash
        for path, digest in plan['pins'].items():
            assert sha(path) == digest, path
        for name, digest in receipt['artifacts'].items():
            assert sha(root / name) == digest, name
    verify_inputs()
    con = sqlite3.connect('file:' + str(Path(plan['bridge']).resolve()) + '?mode=ro', uri=True)
    lookup = {}
    for taxon, protein, model, version in con.execute('SELECT taxon_id,protein_id,model_id,version FROM structures'):
        gene = taxon + '_' + protein
        assert gene not in lookup
        lookup[gene] = (model, version)
    con.close()
    assert len(lookup) == receipt['frozen_model_links']
    checked = []
    for entry in plan['guides']:
        reports = defaultdict(set)
        with Path(entry['events']).open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                reports[row['Orthogroup']].add(row['Gene Tree Node'])
        taxa = {line.split(': ', 1)[1].rsplit('.', 1)[0] for line in Path(entry['species']).read_text().splitlines()}
        counts = Counter(); seen_families = set(); total = 0
        with (root / (entry['guide'] + '_terminal_sisters.tsv')).open() as output, Path(entry['trees']).open() as trees:
            records = iter(csv.DictReader(output, delimiter='\t'))
            current = next(records, None)
            for line in trees:
                family, newick = line.rstrip().split(': ', 1)
                assert family not in seen_families
                seen_families.add(family)
                expected = reconstruct(newick, reports[family], lookup, taxa)
                observed_nodes = set()
                while current is not None and current['family'] == family:
                    assert current['guide'] == entry['guide']
                    name = current['gene_node']
                    assert name in expected and name not in observed_nodes
                    observed_nodes.add(name)
                    compare({k: v for k, v in current.items() if k not in ['guide', 'family']}, expected[name])
                    counts[current['candidate_status'] + ':' + current['model_coverage']] += 1
                    total += 1
                    current = next(records, None)
                assert observed_nodes == set(expected), family
                if len(seen_families) % 1000 == 0:
                    print(entry['guide'], len(seen_families), 'trees verified', total, 'pairs', flush=True)
            assert current is None, 'Extra output records'
        summary = dict(guide=entry['guide'], trees=len(seen_families), pairs=total, counts=dict(counts))
        assert len(seen_families) == entry['expected_trees']
        assert summary == next(row for row in receipt['guides'] if row['guide'] == entry['guide'])
        checked.append(summary)
    verify_inputs(); verify_config()
    result = dict(status='passed_full_terminal_sister_inventory_readback',
                  producer_receipt_sha256=receipt_hash, config_sha256=config_hash,
                  plan_sha256=sha(config['plan']), guides=checked, checker_sha256=sha(__file__),
                  scope='All native trees traversed independently in postorder; exact terminal pair sets, every exported field, model joins, branch lengths and all summary counts reconstructed. No inference of speciation, biological orthology, matched comparability or duplication effect.')
    target = Path(config['output']); assert not target.exists()
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    csv.field_size_limit(32 * 1024 * 1024)
    main()
