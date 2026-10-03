#!/usr/bin/env python3
"""Diagnose the first complete cohort of a finished, unclosed covariance producer.

This is a bounded diagnostic, not full covariance qualification or a refit.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); assert not args.output.exists()
    plan_path = Path('metadata/full_uniform_covariance_qualification_plan_20261002.json')
    plan = json.loads(plan_path.read_text()); root = Path(plan['output'])
    receipt_path = root / 'receipt.json'; receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'complete_full_uniform_covariance_qualification_pending_independent_readback'
    assert receipt['plan_sha256'] == sha(plan_path)
    path = root / 'design_covariance_audits.jsonl.gz'
    assert sha(path) == receipt['artifacts'][path.name]
    rows = []
    with gzip.open(path, 'rt') as handle:
        for _ in range(300): rows.append(json.loads(next(handle)))
    assert len({row['cohort_id'] for row in rows}) == 1
    assert len({row['design_id'] for row in rows}) == 30
    assert Counter(row['loading_mode'] for row in rows) == {'signed': 150, 'unsigned': 150}
    names = ['residual', 'background_node', 'model_pair', 'gene', 'model', 'family_intercept', 'species']
    relations = {
        'pair_minus_residual_minus_background': [-1, -1, 1, 0, 0, 0, 0],
        'gene_minus_half_pair': [0, 0, -.5, 1, 0, 0, 0],
        'model_minus_gene': [0, 0, 0, -1, 1, 0, 0]}
    diagnostics = []
    for row in rows:
        audit = row['numerical_audit']; assert audit['kernel_names'] == names
        record = dict(audit_id=row['audit_id'], design_id=row['design_id'], tree=row['tree'], loading_mode=row['loading_mode'])
        for kind, field, bound_field in [('raw', 'raw_gram', 'raw_roundoff_envelope'),
                                        ('projected', 'projected_gram', 'projected_roundoff_envelope')]:
            gram = np.asarray(audit[field]); envelope = np.asarray(audit[bound_field])
            assert gram.shape == envelope.shape == (7, 7)
            checks = {}
            for name, coefficients in relations.items():
                vector = np.asarray(coefficients)
                error = np.abs(gram @ vector)
                bound = envelope @ np.abs(vector) + 100 * np.finfo(float).eps * np.max(np.abs(gram))
                checks[name] = dict(within_recorded_roundoff=bool(np.all(error <= bound)), max_absolute_gram_residual=float(np.max(error)))
            selected = [0, 1, 5, 6]
            submatrix = gram[np.ix_(selected, selected)]
            diagonal = np.diag(submatrix)
            assert np.all(diagonal > 0)
            normalized = submatrix / np.sqrt(np.outer(diagonal, diagonal))
            record[kind] = dict(relations=checks, selected_four_normalized_gram_singular_values=np.linalg.svd(normalized, compute_uv=False).tolist())
        diagnostics.append(record)
    result = dict(status='completed_first_cohort_covariance_dependency_diagnostic',
        checked_utc=datetime.now(timezone.utc).isoformat(), producer_audit_rows=receipt['audit_rows'],
        producer_status_counts=receipt['audit_status_counts'], producer_setting_links=receipt['setting_audit_links'],
        checked_cohorts=1, checked_designs=30, checked_audits=len(rows), cohort_id=rows[0]['cohort_id'],
        relations=relations, diagnostics=diagnostics, independent_full_readback_complete=False,
        proposed_covariance_model_accepted=False, scientific_eligibility=False,
        source_hashes={str(p): sha(p) for p in [Path(__file__), plan_path, receipt_path, path]},
        scope='Producer summary reports all1,302,000bases require review; every audit in only the first complete cohort is diagnosed here. Gram-vector residuals are compared with recorded roundoff envelopes, not used to assert exact full-grid operator identities. Four-term singular values are descriptive diagnostics, not full qualification. Every original record/job is preserved. A reduced parametrization needs exact cohort/operator and nonnegative variance-cone proofs, full design/REML qualification and calibration before fitting.')
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ['source_hashes', 'diagnostics']}, indent=2))


if __name__ == '__main__':
    main()
