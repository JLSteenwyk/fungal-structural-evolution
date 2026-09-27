#!/usr/bin/env python3
"""Rebuild all ecological bootstrap summary fields with independent pandas joins."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from Bio import Phylo
import numpy as np
import pandas as pd
import psutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--wait-for', type=Path)
    args = parser.parse_args()
    script_hash = sha(__file__)
    if args.wait_for:
        identity = json.loads(args.wait_for.read_text())
        while True:
            try:
                process = psutil.Process(identity['pid'])
                if process.create_time() != identity['created'] or process.status() == psutil.STATUS_ZOMBIE:
                    break
                if process.cmdline() != identity['cmdline']:
                    raise ValueError('Summary controller identity changed')
            except psutil.NoSuchProcess:
                break
            print('waiting_for_complete_ecology_summary', identity['pid'], flush=True)
            time.sleep(30)
    receipt_path = args.input / 'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'complete_audited_ecology_bootstrap_uncertainty_summary_pending_readback'
    for path, digest in receipt['sources'].items():
        assert sha(path) == digest, path
    for name, digest in receipt['artifacts'].items():
        assert sha(args.input / name) == digest, name
    boot = Path('results/ecology/bootstrap-optimal-edge-states-20260927-v1')
    mlroot = Path('results/ecology/optimal-edge-states-20260927-v1')
    read = lambda p: pd.read_csv(p, sep='\t', keep_default_na=False)
    frequencies = read(boot / 'split_frequencies.tsv')
    ml = read(mlroot / 'edge_states.tsv')
    actual = read(args.input / 'split_uncertainty.tsv')
    keys = ['tree_source', 'coding', 'split_id']
    for frame in [frequencies, ml, actual]:
        assert not frame.duplicated(keys).any()
    assert set(map(tuple, actual[keys].values)) == set(map(tuple, frequencies[keys].values)) | set(map(tuple, ml[keys].values))
    source_plan = json.loads(Path('metadata/ecology_bootstrap_edges_plan_20260927.json').read_text())
    taxa = {r['label']: {t.name for t in Phylo.read(r['tree'], 'newick').get_terminals()}
            for r in source_plan['trees']}
    ml_sides = []
    for r in ml.itertuples(index=False):
        side = sorted(r.child_taxa.split(';'))
        other = sorted(taxa[r.tree_source] - set(side))
        canonical = side if (len(side), side) <= (len(other), other) else other
        assert hashlib.sha256(json.dumps(canonical, separators=(',', ':')).encode()).hexdigest() == r.split_id
        ml_sides.append(';'.join(canonical))
    ml = ml.assign(ml_canonical_side=ml_sides)
    expected = frequencies.merge(ml[keys + ['status', 'ml_canonical_side']], on=keys, how='outer', validate='one_to_one')
    both = expected.canonical_side.notna() & expected.ml_canonical_side.notna()
    assert (expected.loc[both, 'canonical_side'] == expected.loc[both, 'ml_canonical_side']).all()
    expected['canonical_side'] = expected.canonical_side.fillna(expected.ml_canonical_side)
    expected['ml_status'] = expected.status.fillna('absent_from_ml_tree')
    countcols = ['present', 'required_change', 'optional_change', 'no_change_in_any_optimum']
    for col in countcols:
        expected[col] = expected[col].fillna(0).astype(int)
    expected['ensemble_trees'] = 1000
    expected['absent_trees'] = 1000 - expected.present
    expected['presence_fraction'] = expected.present / 1000
    for col in countcols[1:]:
        expected[col + '_fraction_all_trees'] = expected[col] / 1000
        expected[col + '_fraction_when_present'] = expected[col] / expected.present.replace(0, np.nan)
    expected = expected.set_index(keys).sort_index()
    actual = actual.set_index(keys).sort_index()
    assert actual.index.equals(expected.index)
    assert set(actual.columns) == set(expected.columns) - {'status', 'ml_canonical_side'}
    for col in actual:
        if col in ['canonical_side', 'ml_status']:
            assert actual[col].equals(expected[col]), col
        else:
            left = pd.to_numeric(actual[col].replace('', np.nan)).to_numpy(float)
            right = expected[col].to_numpy(float)
            np.testing.assert_allclose(left, right, rtol=0, atol=1e-14, equal_nan=True)
    trees = read(boot / 'tree_summaries.tsv')
    metriccols = ['minimum_changes', 'required_change', 'optional_change', 'no_change_in_any_optimum']
    long = trees.melt(id_vars=['tree_source', 'coding', 'tree_index'], value_vars=metriccols, var_name='metric', value_name='value')
    distkeys = ['tree_source', 'coding', 'metric', 'value']
    distributions = long.groupby(distkeys).size().rename('trees').to_frame()
    distributions['ensemble_trees'] = 1000
    distributions['fraction'] = distributions.trees / 1000
    observed = read(args.input / 'tree_metric_distributions.tsv')
    assert not observed.duplicated(distkeys).any()
    pd.testing.assert_frame_equal(observed.set_index(distkeys).sort_index(), distributions.sort_index(), check_dtype=False, atol=1e-14, rtol=0)
    assert len(actual) == receipt['split_rows'] and len(observed) == receipt['distribution_rows']
    assert len(trees) == receipt['coding_tree_combinations'] == 4000
    assert sha(__file__) == script_hash
    result = dict(status='passed_full_ecology_bootstrap_summary_readback',
                  source_receipt_sha256=sha(receipt_path), script_sha256=script_hash,
                  split_rows=len(actual), distribution_rows=len(observed), coding_tree_combinations=len(trees),
                  scope='Independent full pandas outer joins reconstruct counts, both denominator choices, missing conditional fractions, canonical ML split annotations and every per-tree metric distribution. Source network-flow audit establishes underlying edge costs; this check establishes summary aggregation, not biological origins or effects.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
