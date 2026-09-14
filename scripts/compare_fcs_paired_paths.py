#!/usr/bin/env python3
"""Compare baseline and FCS refits on exactly the same accepted taxon pairs."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table

PATHS = ['aa_tree_path_point', '3di_af_tree_path_point', '3di_af_empirical_tree_path_point', '3di_llm_tree_path_point']
METRICS = ['ca_superposition_rmsd_angstrom', 'local_distance_mean_absolute_change_angstrom', 'pae10_local_mean_absolute_change_angstrom']
KEYS = ['marker', 'taxon_a', 'taxon_b']


def rank_correlation(x, y):
    if len(x) < 3 or len(set(x)) < 2 or len(set(y)) < 2:
        return None
    return float(np.corrcoef(pd.Series(x).rank().to_numpy(), pd.Series(y).rank().to_numpy())[0, 1])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['baseline', 'sensitivity', 'inputs', 'output']:
        p.add_argument('--' + key, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    receipts = {k: checked_receipt(getattr(a, k)) for k in ['baseline', 'sensitivity', 'inputs']}
    if any(receipts[k]['status'] != 'complete_paired_site_tree_path_point_benchmark' for k in ['baseline', 'sensitivity']):
        raise ValueError('Incomplete benchmarks')
    if receipts['sensitivity']['source_receipts']['inputs'] != sha(a.inputs / 'receipt.json') or receipts['baseline']['source_receipts']['inputs'] != receipts['inputs']['source_input_receipt_sha256']:
        raise ValueError('Input lineage differs')
    frames = {k: pd.read_csv(getattr(a, k) / 'path_geometry_points.tsv', sep='\t', dtype=str, keep_default_na=False) for k in ['baseline', 'sensitivity']}
    before, after = frames['baseline'], frames['sensitivity']
    if list(before.columns) != list(after.columns) or any(f.duplicated(KEYS).any() for f in frames.values()):
        raise ValueError('Duplicate keys or changed schemas')
    omitted = {(r['marker'], r['taxon_id']) for r in read_table(a.inputs / 'omitted_marker_observations.tsv')}
    markers = set(after.marker)
    mask = [m in markers and (m, x) not in omitted and (m, y) not in omitted for m,x,y in before[KEYS].itertuples(index=False, name=None)]
    before = before.loc[mask].sort_values(KEYS).reset_index(drop=True)
    after = after.sort_values(KEYS).reset_index(drop=True)
    fixed = [c for c in before.columns if c not in PATHS + ['tree_taxa']]
    if not before[fixed].equals(after[fixed]) or len(after) != receipts['sensitivity']['accepted_pairs']:
        raise ValueError('Common pair grid or fixed geometry fields changed')
    disposition = {r['marker']: r for r in read_table(a.inputs / 'marker_disposition.tsv')}
    if len(markers) != receipts['inputs']['ready_markers']:
        raise ValueError('Marker universe differs')
    rows, ranks = [], []
    for marker, group in after.groupby('marker', sort=True):
        original = before.loc[group.index]
        for frame, field in [(original, 'baseline_eligible_taxa'), (group, 'sensitivity_eligible_taxa')]:
            if set(frame.tree_taxa) != {disposition[marker][field]}:
                raise ValueError('Tree taxon counts differ')
        for column in PATHS:
            x, y = original[column].to_numpy(float), group[column].to_numpy(float)
            if not np.isfinite([x,y]).all() or min(x.min(), y.min()) < 0:
                raise ValueError('Invalid path values')
            delta = y-x
            rows.append({'marker': marker, 'path_model': column, 'common_pairs': len(x),
                         'baseline_path_median': float(np.median(x)), 'sensitivity_path_median': float(np.median(y)),
                         'median_signed_path_change': float(np.median(delta)),
                         'median_absolute_path_change': float(np.median(abs(delta))),
                         'maximum_absolute_path_change': float(max(abs(delta))),
                         'baseline_sensitivity_rank_correlation': rank_correlation(x,y)})
            for metric in METRICS:
                geometry = np.array([float(v) if v else np.nan for v in group[metric]])
                valid = np.isfinite(geometry)
                old, new = rank_correlation(x[valid], geometry[valid]), rank_correlation(y[valid], geometry[valid])
                ranks.append({'marker': marker, 'path_model': column, 'geometry_metric': metric,
                              'common_pairs': int(valid.sum()), 'baseline_rank_correlation': old,
                              'sensitivity_rank_correlation': new,
                              'rank_correlation_change': new-old if new is not None and old is not None else None})
    a.output.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(a.output / 'marker_path_changes.tsv', sep='\t', index=False)
    pd.DataFrame(ranks).to_csv(a.output / 'matched_geometry_rank_changes.tsv', sep='\t', index=False)
    result = {'status': 'complete_common_pair_fcs_path_comparison', 'markers': len(markers),
              'common_accepted_pairs': len(after), 'path_values_compared': len(after)*len(PATHS),
              'unchanged_fields_per_pair': len(fixed), 'marker_model_comparisons': len(rows),
              'marker_model_geometry_comparisons': len(ranks),
              'source_receipts': {k: sha(getattr(a,k)/'receipt.json') for k in receipts},
              'script_sha256': sha(Path(__file__)),
              'artifacts': {f.name: sha(f) for f in a.output.iterdir()},
              'interpretation': 'Within-model fitted path changes and geometry rank changes on identical retained pairs. Baseline trees retain flagged observations during fitting; sensitivity trees refit after omission. Geometry and masks are identical for compared pairs. Descriptive, without p-values or independent-pair assumptions. No path ratios, physical displacement interpretation, acceleration test or coupling robustness claim.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
