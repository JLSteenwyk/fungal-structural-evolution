#!/usr/bin/env python3
"""Measure full-family reconciliation dimensions without enumerating protein pairs."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import time
from orthofinder.tools import mcl


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    a = ap.parse_args()
    plan = json.loads(a.plan.read_text())
    for p, expected in plan['pins'].items():
        if sha(Path(p)) != expected:
            raise ValueError('Changed pinned input ' + p)
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    start = time.time()
    universe = json.loads(Path(plan['universe_receipt']).read_text())
    inventory = json.loads(Path(plan['inventory_receipt']).read_text())
    inv_table = Path(plan['inventory_receipt']).parent / 'families.tsv'
    if sha(inv_table) != inventory['artifacts']['families.tsv']:
        raise ValueError('Changed family inventory')
    with inv_table.open() as f:
        sizes = {x['family']: int(x['source_sequences']) for x in csv.DictReader(f, delimiter='\t')}
    clusters = Path(plan['clusters'])
    if sha(clusters) != universe['cluster_sha256']:
        raise ValueError('Cluster source differs from full-universe audit')
    groups = mcl.GetPredictedOGs(str(clusters))
    if len(groups) != universe['counts']['families'] or len(sizes) != len(groups):
        raise ValueError('Family count mismatch')
    records = []
    all_species = set()
    for i, genes in enumerate(groups):
        counts = Counter(g.split('_', 1)[0] for g in genes)
        all_species.update(counts)
        n, k = len(genes), len(counts)
        family = f'OG{i:07d}'
        if sizes.get(family) != n:
            raise ValueError('Native family ordering/count differs from audited input')
        # Two distinct O(number-of-species) calculations avoid materializing pairs.
        pairs = (n * n - sum(v * v for v in counts.values())) // 2
        seen = independent = 0
        for value in counts.values():
            independent += seen * value
            seen += value
        if pairs != independent:
            raise ValueError('Candidate pair arithmetic differs')
        records.append(dict(family=family, genes=n, taxa=k, maximum_copies_per_taxon=max(counts.values()),
                            unordered_cross_species_candidate_pairs=pairs,
                            occupied_unordered_taxon_pairs=k * (k - 1) // 2))
    if sum(x['genes'] for x in records) != universe['counts']['genes'] or len(all_species) != 526:
        raise ValueError('Full scope differs')
    table = out / 'family_workload.tsv'
    with table.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(records)
    tree_dir = Path(plan['tree_directory'])
    trees = sorted(tree_dir.glob('OG*.txt'))
    if {p.stem for p in trees} != {x['family'] for x in records if x['genes'] >= 3}:
        raise ValueError('Expected tree filename universe differs')
    tree_sizes = [(p.name, p.stat().st_size) for p in trees]
    top = sorted(records, key=lambda x: (-x['unordered_cross_species_candidate_pairs'], x['family']))[:25]
    total_pairs = sum(x['unordered_cross_species_candidate_pairs'] for x in records)
    result = dict(status='complete_full_reconciliation_workload_inventory',
                  families=len(records), genes=sum(x['genes'] for x in records), taxa=len(all_species),
                  tree_eligible_families=sum(x['genes'] >= 3 for x in records),
                  single_taxon_families=sum(x['taxa'] == 1 for x in records),
                  unordered_cross_species_candidate_pairs=total_pairs,
                  directed_cross_species_candidate_pairs=2 * total_pairs,
                  occupied_unordered_family_taxon_pairs=sum(x['occupied_unordered_taxon_pairs'] for x in records),
                  theoretical_directed_species_pair_files=526 * 525,
                  top25_families=top,
                  top25_candidate_pair_fraction=sum(x['unordered_cross_species_candidate_pairs'] for x in top) / total_pairs,
                  source_tree_files=len(trees), source_tree_bytes_at_observation=sum(n for _, n in tree_sizes),
                  empty_source_trees_at_observation=[name for name, n in tree_sizes if n == 0],
                  flat_directed_pair_storage_scenarios=[dict(bytes_per_pair=b, gib=2 * total_pairs * b / 2**30) for b in [64, 128, 256]],
                  plan_sha256=sha(a.plan), script_sha256=sha(Path(__file__)), elapsed_seconds=time.time() - start,
                  artifacts={'family_workload.tsv': sha(table)},
                  interpretation='Exact candidate comparison dimensions from audited family membership; pairs are not inferred orthologs and are not materialized. Flat-table storage scenarios are arithmetic illustrations, not predictions of grouped native output. Live tree byte counts are an observation excluding unfinished repair content. These dimensions do not establish runtime or peak memory; final launch resource plan remains required.')
    for p, expected in plan['pins'].items():
        if sha(Path(p)) != expected:
            raise ValueError('Pinned input changed during inventory')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'top25_families'}, indent=2))


if __name__ == '__main__':
    main()
