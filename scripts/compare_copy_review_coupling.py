#!/usr/bin/env python3
"""Verify and compare full-cohort and explicit copy-review omission coupling fits."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def checked(folder):
    receipt = json.loads((folder / 'receipt.json').read_text())
    for name, digest in receipt['artifacts'].items():
        if sha(folder / name) != digest:
            raise ValueError('Changed artifact ' + str(folder / name))
    return receipt


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['full', 'omission', 'full-resampling', 'omission-resampling', 'output']:
        ap.add_argument('--' + name, type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    full, omission = checked(a.full), checked(a.omission)
    rf, ro = checked(a.full_resampling), checked(a.omission_resampling)
    assert rf['source_fit_receipt_sha256'] == sha(a.full / 'receipt.json')
    assert ro['source_fit_receipt_sha256'] == sha(a.omission / 'receipt.json')
    assert full['source_frame_sha256'] == omission['source_frame_sha256']
    assert full['source_accessibility_sha256'] == omission['source_accessibility_sha256']
    pf, po = full['marker_review_policy'], omission['marker_review_policy']
    assert pf['mode'] == 'retain_flagged_exploratory' and po['mode'] == 'exclude_reviewed_markers'
    assert pf['expected_flagged_markers'] == po['expected_flagged_markers']
    assert len(pf['expected_flagged_markers']) == 1
    marker = pf['expected_flagged_markers'][0]['marker']
    data = pd.read_csv(a.full / 'site_covariates.tsv', sep='\t')
    reduced = pd.read_csv(a.omission / 'site_covariates.tsv', sep='\t')
    expected = data[data.marker != marker].reset_index(drop=True)
    pd.testing.assert_frame_equal(expected, reduced, check_exact=True)
    assert len(data) - len(reduced) == pf['expected_flagged_markers'][0]['sites']
    assert reduced.marker.nunique() == data.marker.nunique() - 1
    cf = pd.read_csv(a.full / 'focal_coefficients.tsv', sep='\t')
    co = pd.read_csv(a.omission / 'focal_coefficients.tsv', sep='\t')
    assert len(cf) == len(co) == 72
    matched = cf.merge(co, on=['model_id', 'term'], suffixes=('_full', '_omission'), validate='one_to_one')
    assert len(matched) == 72
    influence = pd.read_csv(a.full_resampling / 'marker_influence.tsv', sep='\t')
    aliases = {'aa_log1p_rate': 'aa_log1p_rate_at_RSA_0.25', 'rsa_centered': 'rsa_at_aa_log1p_rate_0', 'aa_by_rsa': 'interaction'}
    check = co[['model_id', 'term', 'coefficient']].copy()
    check['contrast'] = check.term.map(aliases)
    check = check.merge(influence[influence.omitted_marker == marker], on=['model_id', 'contrast'], validate='one_to_one')
    assert len(check) == 72
    assert np.allclose(check.coefficient, check.estimate_without_marker, rtol=1e-7, atol=1e-9)
    matched['omission_minus_full'] = matched.coefficient_omission - matched.coefficient_full
    summary = []
    bootstrap = []
    for label, coefficients, folder, receipt in [('full', cf, a.full_resampling, rf), ('omission', co, a.omission_resampling, ro)]:
        assert receipt['status'] == 'complete_conditional_marker_resampling'
        assert receipt['models'] == 24 and receipt['bootstrap_fits'] == 48000
        assert receipt['leave_one_marker_out_fits'] == 24 * receipt['markers']
        for term, group in coefficients.groupby('term'):
            summary.append(dict(analysis=label, term=term, specifications=len(group),
                                positive_coefficients=int((group.coefficient > 0).sum()),
                                bh_q_below_005=int((group.bh_q_across_all_72_focal_tests < .05).sum())))
        contrasts = pd.read_csv(folder / 'contrast_summary.tsv', sep='\t')
        assert len(contrasts) == 120
        for contrast, group in contrasts.groupby('contrast'):
            bootstrap.append(dict(analysis=label, contrast=contrast, specifications=len(group),
                                  intervals_above_zero=int((group.bootstrap_percentile_lower > 0).sum()),
                                  intervals_below_zero=int((group.bootstrap_percentile_upper < 0).sum())))
    a.output.mkdir(parents=True)
    matched.to_csv(a.output / 'matched_focal_coefficients.tsv', sep='\t', index=False)
    pd.DataFrame(summary).to_csv(a.output / 'term_summary.tsv', sep='\t', index=False)
    pd.DataFrame(bootstrap).to_csv(a.output / 'bootstrap_summary.tsv', sep='\t', index=False)
    result = dict(status='passed_copy_review_coupling_comparison',
                  source_receipts={str(p): sha(p / 'receipt.json') for p in [a.full, a.omission, a.full_resampling, a.omission_resampling]},
                  omitted_marker=marker, full_markers=full['markers'], omission_markers=omission['markers'],
                  full_sites=len(data), omission_sites=len(reduced), matched_focal_coefficients=72,
                  full_loo_vs_separate_omission_max_difference=float(np.max(np.abs(check.coefficient - check.estimate_without_marker))),
                  script_sha256=sha(Path(__file__)), artifacts={p.name: sha(p) for p in a.output.iterdir()},
                  interpretation='Exact covariate subset and all 72 focal omission coefficients checked against independently absorbed full-cohort leave-one-marker-out estimates. BH families remain separate per analysis; bootstrap intervals are unadjusted sensitivity intervals. Results remain conditional on estimated rates, fixed gene trees and predicted structures; copy review is unresolved and no causal or confirmatory single-ortholog conclusion follows.')
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
