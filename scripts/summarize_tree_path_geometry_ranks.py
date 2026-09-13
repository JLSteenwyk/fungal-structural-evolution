#!/usr/bin/env python3
"""Describe within-marker path/geometry rank associations without significance tests."""
import argparse
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def rank_association(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    return float(np.corrcoef(rankdata(x, method='average'), rankdata(y, method='average'))[0, 1])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['benchmark', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use immutable output')
    r = checked_receipt(a.benchmark)
    groups = defaultdict(list)
    for row in read_table(a.benchmark / 'path_geometry_points.tsv'):
        groups[row['marker']].append(row)
    if sum(map(len, groups.values())) != r['accepted_pairs']:
        raise ValueError('Input row count mismatch')
    geometry = ['ca_superposition_rmsd_angstrom', 'local_distance_mean_absolute_change_angstrom', 'pae10_local_mean_absolute_change_angstrom']
    paths = ['aa_tree_path_point', '3di_af_tree_path_point', '3di_af_empirical_tree_path_point', '3di_llm_tree_path_point']
    results = []
    for marker, rows in sorted(groups.items()):
        for geom in geometry:
            cohort = [row for row in rows if row[geom] != '']
            for path in paths:
                x, y = [float(row[path]) for row in cohort], [float(row[geom]) for row in cohort]
                if not all(np.isfinite(x)) or not all(np.isfinite(y)):
                    raise ValueError('Nonfinite input')
                rho = rank_association(x, y)
                results.append({'marker': marker, 'path_metric': path, 'geometry_metric': geom,
                                'pairs': len(cohort), 'omitted_geometry_pairs': len(rows)-len(cohort),
                                'rank_correlation': '' if rho is None else rho,
                                'status': 'descriptive_estimable' if rho is not None else 'insufficient_pairs_or_constant_values',
                                'marker_review_status': rows[0]['marker_review_status']})
    summaries = []
    for path in paths:
        for geom in geometry:
            selected = [row for row in results if row['path_metric'] == path and row['geometry_metric'] == geom and row['rank_correlation'] != '']
            values = [row['rank_correlation'] for row in selected]
            summaries.append({'path_metric': path, 'geometry_metric': geom, 'estimable_markers': len(values),
                              'median_marker_rank_correlation': float(np.median(values)) if values else '',
                              'negative_marker_correlations': sum(v < 0 for v in values)})
    a.output.mkdir(parents=True)
    write_table(a.output / 'marker_rank_associations.tsv', results)
    write_table(a.output / 'equal_marker_summary.tsv', summaries)
    out = {'status': 'complete_descriptive_within_marker_rank_summary', 'source_receipt_sha256': sha(a.benchmark / 'receipt.json'),
           'script_sha256': sha(Path(__file__)), 'markers': len(groups), 'summary_rows': len(results),
           'interpretation': 'Average-tie rank correlations within markers, no p-values or independence claim. Equal-marker summaries are descriptive, not meta-analysis. Geometry availability defines each metric cohort; all four models use that same cohort. Shared ancestry, topology/branch uncertainty, source/confidence selection and model adequacy remain unadjusted. Correlation does not establish additive physical distances or evolutionary sequence–structure coupling.',
           'artifacts': {f.name: sha(f) for f in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(out, indent=2) + '\n'); print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
