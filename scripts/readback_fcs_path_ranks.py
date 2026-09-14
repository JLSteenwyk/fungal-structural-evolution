#!/usr/bin/env python3
"""Recalculate all reported common-pair rank correlations using SciPy."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['baseline', 'sensitivity', 'comparison', 'output']:
        p.add_argument('--' + key, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    r = checked_receipt(a.comparison)
    keys = ['marker', 'taxon_a', 'taxon_b']
    frames = {}
    for key in ['baseline', 'sensitivity']:
        folder = getattr(a, key); checked_receipt(folder)
        if r['source_receipts'][key] != sha(folder / 'receipt.json'):
            raise ValueError('Source differs')
        # Preserve IEEE floats: default CSV parsing can create/break rank ties.
        frame = pd.read_csv(folder / 'path_geometry_points.tsv', sep='\t', float_precision='round_trip').set_index(keys)
        if not frame.index.is_unique:
            raise ValueError('Duplicate source keys')
        frames[key] = frame
    new = frames['sensitivity']; old = frames['baseline'].loc[new.index]
    differences = []
    def check(x, y, reported):
        observed = float(spearmanr(x, y).statistic)
        if np.isnan(observed) and pd.isna(reported):
            return
        error = abs(observed - reported)
        if not np.isfinite(error) or error > 1e-12:
            raise ValueError('Independent rank correlation differs')
        differences.append(error)
    for row in pd.read_csv(a.comparison / 'marker_path_changes.tsv', sep='\t', float_precision='round_trip').to_dict('records'):
        marker, column = row['marker'], row['path_model']
        check(old.loc[marker, column], new.loc[marker, column], row['baseline_sensitivity_rank_correlation'])
    for row in pd.read_csv(a.comparison / 'matched_geometry_rank_changes.tsv', sep='\t', float_precision='round_trip').to_dict('records'):
        marker, column, metric = row['marker'], row['path_model'], row['geometry_metric']
        geometry = new.loc[marker, metric]
        valid = np.isfinite(geometry)
        for frame, label in [(old, 'baseline'), (new, 'sensitivity')]:
            check(frame.loc[marker, column][valid], geometry[valid], row[label + '_rank_correlation'])
        delta = row['sensitivity_rank_correlation'] - row['baseline_rank_correlation']
        if not np.isclose(delta, row['rank_correlation_change'], atol=1e-12, rtol=0, equal_nan=True):
            raise ValueError('Rank change differs')
    result = {'status': 'passed_full_reported_rank_correlation_readback',
              'common_pairs': len(new), 'correlations_checked': len(differences),
              'maximum_absolute_difference': max(differences),
              'comparison_receipt_sha256': sha(a.comparison / 'receipt.json'),
              'script_sha256': sha(Path(__file__)),
              'scope': 'SciPy Spearman recalculation of all path-ranking and before/after geometry-ranking correlations, joined independently by taxon-pair identity. No p-values used. Does not independently refit trees or estimate uncertainty.'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
