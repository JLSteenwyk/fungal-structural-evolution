#!/usr/bin/env python3
"""Build domain-averaged polynomial identity contrasts, checking SQL with pandas."""
import argparse
import json
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    plan_hash = sha(args.plan)

    def verify():
        assert sha(args.plan) == plan_hash
        for path, digest in plan['pins'].items():
            assert sha(path) == digest, path

    verify()
    source = Path(plan['summaries'])
    proof = json.loads(Path(plan['summary_audit']).read_text())
    assert proof['status'] == 'passed_full_matched_domain_record_summary_readback'
    assert proof['source_receipt_sha256'] == sha(source / 'receipt.json')
    measurement_proof = json.loads(Path(plan['measurement_audit']).read_text())
    assert measurement_proof['source_receipt_sha256'] == sha(Path(plan['measurements']).parent / 'receipt.json')
    assert measurement_proof['status'] == plan['measurement_audit_status']
    parts = json.loads((source / 'partition_manifest.json').read_text())
    assert len(parts) == 192
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    db = duckdb.connect()
    db.execute("SET threads=1")
    db.execute("SET memory_limit='3GB'")
    db.execute("CREATE TABLE raw AS SELECT * FROM read_csv(?, delim='\t', header=true, all_varchar=true)", [plan['measurements']])
    # Independent parser and arithmetic for every output, not a sample.
    raw = pd.read_csv(plan['measurements'], sep='\t', dtype=str, keep_default_na=False)
    assert len(raw) == db.execute('SELECT count(*) FROM raw').fetchone()[0] == 1200640
    checked = 0
    artifacts = {}
    diagnostics = []
    for part in parts:
        original_path = source / part['path']
        assert sha(original_path) == part['sha256']
        flag = part['screen'] + ('_mask_pass' if part['cohort'] == 'same_mask' else '_both_masks_pass')
        keys = ['boundary', 'mask', 'target_order', 'background_order']
        where = ' AND '.join(k + '=?' for k in keys) + f" AND {flag}='1'"
        x = 'CAST(target_sequence_identity_exact AS DOUBLE)'
        y = 'CAST(background_sequence_identity_exact AS DOUBLE)'
        expressions = [f'avg(pow({x},{p})-pow({y},{p})) AS identity_power_{p}_difference' for p in [1, 2, 3]]
        frame = db.execute('SELECT domain_config_id,count(*) eligible_domains,' + ','.join(expressions) + ' FROM raw WHERE ' + where + ' GROUP BY domain_config_id ORDER BY domain_config_id', [part[k] for k in keys]).df()
        keep = raw[flag].eq('1')
        for k in keys:
            keep &= raw[k].eq(part[k])
        selected = raw.loc[keep]
        a = selected.target_sequence_identity_exact.astype(float).to_numpy()
        b = selected.background_sequence_identity_exact.astype(float).to_numpy()
        assert np.isfinite(a).all() and np.isfinite(b).all()
        assert ((a >= 0) & (a <= 1) & (b >= 0) & (b <= 1)).all()
        independent = pd.DataFrame({'domain_config_id': selected.domain_config_id.to_numpy(), 'eligible_domains': 1})
        for p in [1, 2, 3]:
            independent[f'identity_power_{p}_difference'] = a ** p - b ** p
        independent = independent.groupby('domain_config_id', sort=True).agg({'eligible_domains': 'sum', **{f'identity_power_{p}_difference': 'mean' for p in [1, 2, 3]}}).reset_index()
        pd.testing.assert_frame_equal(frame, independent, check_dtype=False, atol=1e-12, rtol=1e-12)
        original = pd.read_parquet(original_path).sort_values('domain_config_id').reset_index(drop=True)
        assert frame.domain_config_id.tolist() == original.domain_config_id.tolist()
        np.testing.assert_array_equal(frame.eligible_domains, original.eligible_domains)
        np.testing.assert_allclose(frame.identity_power_1_difference, original.identity_difference, atol=1e-12, rtol=1e-12)
        # Quantify why transforming already averaged identities is insufficient.
        shortcut = original.target_identity ** 2 - original.background_identity ** 2
        error = np.abs(frame.identity_power_2_difference.to_numpy() - shortcut.to_numpy())
        diagnostics.append(dict(index=part['index'], configurations=len(frame), multi_domain_configurations=int((frame.eligible_domains > 1).sum()), quadratic_shortcut_discrepancies_gt_1e_10=int((error > 1e-10).sum()), maximum_quadratic_shortcut_error=float(error.max())))
        path = out / f"{part['index']:03d}.parquet"
        frame.to_parquet(path, index=False)
        pd.testing.assert_frame_equal(pd.read_parquet(path), frame)
        artifacts[path.name] = sha(path)
        checked += len(frame)
        print(f"Verified setting {part['index'] + 1}/192", flush=True)
    db.close()
    verify()
    (out / 'diagnostics.json').write_text(json.dumps(diagnostics, indent=2) + '\n')
    artifacts['diagnostics.json'] = sha(out / 'diagnostics.json')
    receipt = dict(status='complete_polynomial_identity_contrasts_full_dual_implementation_check', plan_sha256=plan_hash, settings=192, configuration_rows=checked, powers=[1, 2, 3], artifacts=artifacts, scope='Mean over eligible domains of identity_target**p minus identity_background**p. All rows checked with separately parsed pandas calculations, serialized readback and original linear contrasts. These are candidate fixed-effect inputs, not fitted nonlinear models or evidence of excess structural change. Model rank, support, convergence and uncertainty require new checks.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
