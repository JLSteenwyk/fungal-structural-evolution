#!/usr/bin/env python3
"""Independently replay control summary medians, quantiles and dispositions with pandas."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_artifacts(folder):
    receipt = json.loads((folder / 'receipt.json').read_text())
    for name, digest in receipt['artifacts'].items():
        if sha(folder / name) != digest:
            raise ValueError(f'Artifact mismatch: {folder / name}')
    return receipt


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['comparisons', 'summary', 'crosswalk', 'output']:
        ap.add_argument('--' + name, type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    cr, sr, xr = [check_artifacts(x) for x in [a.comparisons, a.summary, a.crosswalk]]
    assert sr['comparison_receipt_sha256'] == sha(a.comparisons / 'receipt.json')
    assert cr['source_receipts']['crosswalk'] == sha(a.crosswalk / 'receipt.json')
    def read(folder, name):
        return pd.read_csv(folder / name, sep='\t', dtype={'deposited_model': str})
    data = read(a.comparisons, 'comparisons.tsv')
    rejected = read(a.comparisons, 'exclusions.tsv')
    base = ['rmsd_angstrom', 'local_mean_absolute_difference_angstrom', 'local_rms_difference_angstrom']
    metrics = [p + '_' + m for p in ['af_experiment', 'esm_experiment', 'af_esm'] for m in base]
    deltas = ['paired_delta_' + m for m in base]
    for m in base:
        data['paired_delta_' + m] = data['esm_experiment_' + m] - data['af_experiment_' + m]
    columns = metrics + deltas + ['fraction_full_sequence']
    keys = ['sequence_sha256', 'joint_predicted_plddt_cutoff']
    current = data
    values_checked = 0
    for extra, name in [(['entry_id', 'deposited_model'], 'deposited_model_summary.tsv'), (['entry_id'], 'entry_summary.tsv'), ([], 'protein_summary.tsv')]:
        group = current.groupby(keys + extra, dropna=False)
        expected = group[columns].median().sort_index()
        actual = read(a.summary, name).set_index(keys + extra).sort_index()
        assert expected.index.equals(actual.index), name
        np.testing.assert_allclose(expected.to_numpy(), actual[columns].to_numpy(), rtol=1e-12, atol=1e-12, equal_nan=True)
        np.testing.assert_array_equal(group.size().sort_index().to_numpy(), actual['contributing_units'].to_numpy())
        values_checked += expected.size
        current = expected.reset_index()
    proteins = current
    common = set.intersection(*[set(proteins.loc[proteins[keys[1]] == c, keys[0]]) for c in [0, 70, 90]])
    cohort = read(a.summary, 'cohort_summary.tsv')
    expected_grid = {(label, c, m) for label in ['all_available', 'same_proteins_all_thresholds'] for c in [0,70,90] for m in metrics + deltas}
    assert len(cohort) == len(expected_grid)
    assert set(zip(cohort.cohort, cohort[keys[1]], cohort.metric)) == expected_grid
    quantiles = 0
    for row in cohort.to_dict('records'):
        subset = proteins[proteins[keys[1]] == row[keys[1]]]
        if row['cohort'] == 'same_proteins_all_thresholds':
            subset = subset[subset[keys[0]].isin(common)]
        vals = subset[row['metric']].dropna()
        assert row['proteins'] == len(subset) and row['proteins_with_metric'] == len(vals)
        np.testing.assert_allclose(vals.quantile([.5,.25,.75]).to_numpy(), [row[k] for k in ['median','q25','q75']], rtol=1e-12, atol=1e-12, equal_nan=True)
        quantiles += 3
        if row['metric'] in deltas:
            assert [row[k] for k in ['positive','negative','approximately_zero']] == [sum(vals > 1e-8),sum(vals < -1e-8),sum(abs(vals) <= 1e-8)]
    universe = set(read(a.crosswalk, 'crosswalk.tsv')[keys[0]])
    dispositions = read(a.summary, 'protein_disposition.tsv')
    assert len(dispositions) == len(universe) * 3
    assert set(zip(dispositions[keys[0]], dispositions[keys[1]])) == {(s,c) for s in universe for c in [0,70,90]}
    for row in dispositions.to_dict('records'):
        accepted = data[(data[keys[0]] == row[keys[0]]) & (data[keys[1]] == row[keys[1]])]
        excluded = rejected[(rejected[keys[0]] == row[keys[0]]) & (rejected[keys[1]] == row[keys[1]])]
        assert row['accepted_chain_model_rows'] == len(accepted)
        assert row['excluded_chain_model_rows'] == len(excluded)
        assert row['accepted_entries'] == accepted.entry_id.nunique()
        assert row['status'] == ('included' if len(accepted) else 'no_eligible_chain_model')
    assert sr['control_sequences'] == xr['control_sequences'] == len(universe)
    assert sr['same_proteins_all_thresholds'] == len(common)
    result = {'status':'passed_hierarchical_control_summary_readback', 'hierarchical_metric_values_checked':values_checked, 'cohort_quantiles_checked':quantiles, 'disposition_rows_checked':len(dispositions), 'common_proteins':len(common), 'summary_receipt_sha256':sha(a.summary/'receipt.json'), 'script_sha256':sha(Path(__file__)), 'method':'Independent pandas hierarchical medians, contribution counts, cohort quantiles and delta signs; full crosswalk universe and accepted/excluded disposition counts. Geometry validation remains in the separate coordinate audit.'}
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
