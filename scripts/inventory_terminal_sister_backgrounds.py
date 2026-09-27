#!/usr/bin/env python3
"""Inventory every bifurcating terminal sister pair for background assessment."""
import argparse
import csv
import json
import shutil
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from inventory_duplication_sister_references import tree_index
from run_ortholog_pair_guide_comparison import sha


def pairs(newick, reported, models, taxa):
    nodes, parents, _ = tree_index(newick)
    for name, node in nodes.items():
        if len(node.clades) != 2 or any(child.clades for child in node.clades):
            continue
        a, b = sorted(child.name for child in node.clades)
        ta, tb = a.split('_', 1)[0], b.split('_', 1)[0]
        if ta not in taxa or tb not in taxa:
            raise ValueError('Unknown taxon prefix')
        status = ('reported_duplication' if name in reported else
                  'same_taxon_unreported' if ta == tb else
                  'cross_taxon_unreported_candidate')
        ma, mb = models.get(a), models.get(b)
        coverage = ('identical_model' if ma is not None and ma == mb else
                    'two_distinct_models' if ma is not None and mb is not None else
                    'one_model' if ma is not None or mb is not None else 'no_models')
        yield dict(gene_node=name, gene_a=a, gene_b=b, taxon_a=ta, taxon_b=tb,
                   node_reported_duplication=int(name in reported),
                   parent_node=parents[name] or '',
                   parent_reported_duplication=int(parents[name] in reported),
                   sequence_tip_a=nodes[a].branch_length,
                   sequence_tip_b=nodes[b].branch_length,
                   sequence_pair_distance=nodes[a].branch_length + nodes[b].branch_length,
                   candidate_status=status, model_coverage=coverage,
                   model_a=ma[0] if ma else '', version_a=ma[1] if ma else '',
                   model_b=mb[0] if mb else '', version_b=mb[1] if mb else '')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text()); plan_hash = sha(args.plan)
    def verify():
        assert sha(args.plan) == plan_hash
        for path, digest in plan['pins'].items():
            assert sha(path) == digest, path
    verify()
    out = Path(plan['output'])
    assert shutil.disk_usage(out.parent).free >= plan['resources']['minimum_free_disk_bytes']
    out.mkdir(exist_ok=False)
    con = sqlite3.connect('file:' + str(Path(plan['bridge']).resolve()) + '?mode=ro', uri=True)
    models = {}
    for taxon, protein, model, version in con.execute('SELECT taxon_id,protein_id,model_id,version FROM structures'):
        gene = taxon + '_' + protein
        assert gene not in models
        models[gene] = (model, version)
    con.close()
    summaries = []
    for entry in plan['guides']:
        taxa = {line.split(': ', 1)[1].rsplit('.', 1)[0] for line in Path(entry['species']).read_text().splitlines()}
        assert len(taxa) == 526 and not any('_' in taxon for taxon in taxa)
        reported = defaultdict(set)
        with Path(entry['events']).open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                reported[row['Orthogroup']].add(row['Gene Tree Node'])
        counts = Counter(); families = set(); count = 0; writer = None
        with Path(entry['trees']).open() as handle, (out / (entry['guide'] + '_terminal_sisters.tsv')).open('w') as output:
            for line in handle:
                family, newick = line.rstrip().split(': ', 1)
                assert family not in families
                families.add(family)
                for row in pairs(newick, reported[family], models, taxa):
                    record = dict(guide=entry['guide'], family=family, **row)
                    if writer is None:
                        writer = csv.DictWriter(output, fieldnames=list(record), delimiter='\t', lineterminator='\n')
                        writer.writeheader()
                    writer.writerow(record); count += 1
                    counts[row['candidate_status'] + ':' + row['model_coverage']] += 1
                if len(families) % 1000 == 0:
                    print(entry['guide'], len(families), 'trees', count, 'pairs', flush=True)
        assert len(families) == entry['expected_trees']
        summaries.append(dict(guide=entry['guide'], trees=len(families), pairs=count, counts=dict(counts)))
        print(json.dumps(summaries[-1]), flush=True)
    verify()
    receipt = dict(status='complete_terminal_sister_inventory_pending_independent_readback',
                   plan_sha256=plan_hash, guides=summaries, frozen_model_links=len(models),
                   artifacts={path.name: sha(path) for path in out.glob('*.tsv')},
                   scope='All terminal bifurcating sister pairs in both complete resolved-tree sets. Duplication reports, same-taxon pairs, missing models and identical models retained. Cross-taxon unreported is a candidate label, not verified speciation, orthology or a matched nonduplication control. Guide agreement, orthology membership, matching and numerical structural comparisons remain pending.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == '__main__':
    csv.field_size_limit(32 * 1024 * 1024)
    main()
