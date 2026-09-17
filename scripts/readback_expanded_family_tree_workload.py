#!/usr/bin/env python3
"""Independently verify the task census with the installed native MCL reader."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from orthofinder.tools import mcl


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--merged', type=Path, required=True)
    ap.add_argument('--census', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    receipt = json.loads((a.census / 'receipt.json').read_text())
    if sha(a.merged / 'receipt.json') != receipt['merge_receipt_sha256']:
        raise ValueError('Changed merge receipt')
    for name, digest in receipt['artifacts'].items():
        if sha(a.census / name) != digest:
            raise ValueError('Changed task table')
    with (a.census / 'unique_tree_tasks.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    tasks = {r['membership_sha256']: r for r in rows}
    assert len(tasks) == len(rows) == receipt['unique_tree_tasks']
    observed = Counter()
    merge = json.loads((a.merged / 'receipt.json').read_text())
    for summary in receipt['guides']:
        guide = summary['guide']
        source = next(g for g in merge['guides'] if g['guide'] == guide)
        for name, digest in source['artifacts'].items():
            assert sha(a.merged / name) == digest
        groups = mcl.GetPredictedOGs(str(a.merged / guide / 'clusters_id_pairs.txt'))
        counts = Counter()
        with (a.merged / guide / 'family_sources.tsv').open() as handle:
            crosswalk = list(csv.DictReader(handle, delimiter='\t'))
        assert len(groups) == len(crosswalk) == summary['families']
        for i, genes in enumerate(groups):
            row = crosswalk[i]
            n = len(genes)
            taxa = len({g.split('_', 1)[0] for g in genes})
            bucket = 'singleton' if n == 1 else 'pair' if n == 2 else 'tree'
            counts[row['source_type'] + '_' + bucket + '_families'] += 1
            counts[row['source_type'] + '_' + bucket + '_proteins'] += n
            counts['single_taxon_families'] += taxa == 1
            counts['proteins'] += n
            if n < 3:
                continue
            digest = hashlib.sha256(('\n'.join(sorted(genes)) + '\n').encode()).hexdigest()
            task = tasks[digest]
            assert int(task['proteins']) == n and int(task['taxa']) == taxa
            assert task['source_type'] == row['source_type']
            assert task['original_tree_candidate'] == row['original_tree_candidate']
            if task['representative_guide'] == guide:
                assert task['representative_family'] == f'OG{i:07d}'
                assert task['source_clade'] == row['source_clade']
                assert task['source_family'] == row['source_family']
            observed[digest] += 1
        assert all(summary[k] == v for k, v in counts.items())
        del groups, crosswalk
    assert set(observed) == set(tasks)
    assert all(observed[k] == int(v['guide_occurrences']) for k, v in tasks.items())
    assert sum(v['source_type'] == 'retained' for v in rows) == receipt['retained_tree_candidates']
    assert sum(v['source_type'] == 'discovery' for v in rows) == receipt['unique_new_tree_tasks']
    assert sum(v == 2 for v in observed.values()) == receipt['tasks_shared_across_guides']
    assert max(int(v['proteins']) for v in rows if v['source_type'] == 'discovery') == receipt['largest_new_family']
    result = dict(status='passed_native_expanded_tree_workload_readback',
                  task_rows=len(rows), guide_occurrences=sum(observed.values()),
                  census_receipt_sha256=sha(a.census / 'receipt.json'),
                  script_sha256=sha(Path(__file__)),
                  scope='Every native tree-eligible family represented exactly; sizes, taxa, source crosswalk, guide multiplicity and aggregate counts checked. No tree inference or biological homology validation.')
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
