"""Preserve original estimates and add explicitly selected refined estimates."""
import argparse
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import pandas as pd
import psutil
from ancestral_chain_attempt import sha, write_json
from export_matched_model_grid import COVARIATES

KEY = ['fit_input_id', 'tree']
VALUES = ['intercept', 'conditional_intercept_variance', 'profiled_scale',
          'negative_profiled_reml'] + ['coefficient_'+n for n in COVARIATES] + [
          'variance_ratio_'+n for n in ['background', 'family_component', 'species']]


def overlay(original, refinements):
    assert not original.duplicated(KEY).any()
    result = original.copy().set_index(KEY)
    for name in VALUES:
        result['selected_'+name] = result[name]
    result['selected_source_sha256'] = result.source_fit_sha256
    result['selection'] = 'original_not_targeted'
    result['selected_review_required'] = result.numerical_review_required
    result['refinement_status'] = 'not_targeted'
    result['refinement_source_sha256'] = ''
    seen = set()
    for entry, saved in refinements:
        key = entry['fit_input_id'], entry['tree']
        assert key not in seen and key in result.index
        seen.add(key)
        old = result.loc[key]
        assert bool(old.numerical_review_required)
        assert saved['source']['fit_input_id'] == key[0] and saved['source']['tree'] == key[1]
        assert saved['source']['source_fit_sha256'] == old.source_fit_sha256
        p = saved['payload']
        assert entry['status'] == p['status']
        result.loc[key, 'refinement_status'] = p['status']
        result.loc[key, 'refinement_source_sha256'] = entry['sha256']
        if p['status'] == 'refinement_error_requires_review':
            result.loc[key, 'selection'] = 'original_after_refinement_error'
            result.loc[key, 'selected_review_required'] = True
            continue
        assert p['status'] in ['refinement_passed_numerical_checks', 'refinement_requires_review']
        assert p['records'] == old.records
        active = p['active_covariates']
        assert len(active) == 4 and all(type(v) is bool for v in active)
        assert active == [bool(pd.notna(old['coefficient_'+name])) for name in COVARIATES]
        assert len(p['raw_unit_beta']) == 1 + sum(active)
        cov = np.asarray(p['raw_unit_conditional_beta_covariance'])
        assert cov.shape == (1+sum(active), 1+sum(active)) and np.isfinite(cov).all() and cov[0, 0] > 0
        assert np.isfinite(p['raw_unit_beta']).all()
        assert p['negative_profiled_reml'] <= old.negative_profiled_reml + 1e-7
        values = dict(intercept=p['raw_unit_beta'][0], conditional_intercept_variance=cov[0, 0],
                      profiled_scale=p['profiled_scale'], negative_profiled_reml=p['negative_profiled_reml'])
        beta = iter(p['raw_unit_beta'][1:])
        for name, included in zip(COVARIATES, active):
            values['coefficient_'+name] = next(beta) if included else np.nan
        assert len(p['ratios']) == 3 and np.isfinite(p['ratios']).all()
        for name, value in zip(['background', 'family_component', 'species'], p['ratios']):
            values['variance_ratio_'+name] = value
        for name, value in values.items():
            result.loc[key, 'selected_'+name] = value
        result.loc[key, 'selection'] = 'refined_candidate'
        result.loc[key, 'selected_source_sha256'] = entry['sha256']
        result.loc[key, 'selected_review_required'] = p['status'] != 'refinement_passed_numerical_checks'
    return result.reset_index(), seen


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    config = json.loads(args.plan.read_text()); plan_hash = sha(args.plan)
    def verify():
        assert sha(args.plan) == plan_hash
        for path, digest in config['pins'].items():
            assert sha(path) == digest, path
    verify()
    launch = json.loads(Path(config['audit_launch']).read_text())
    while True:
        state = dict(line.split('=', 1) for line in subprocess.check_output([
            'systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState',
            '-p', 'MainPID', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
        if state['ActiveState'] in ['inactive', 'failed']:
            break
        assert state['ActiveState'] == 'active' and int(state['MainPID']) == launch['pid']
        try:
            process = psutil.Process(launch['pid'])
            assert process.create_time() == launch['created'] and process.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            time.sleep(1)
            continue
        time.sleep(30)
    assert state['ActiveState'] == 'inactive' and state['Result'] == 'success' and state['ExecMainStatus'] == '0'
    verify()
    roots = {name: Path(config[name]) for name in ['original', 'refinement', 'audit']}
    receipts = {}; bindings = {}
    for name, root in roots.items():
        receipts[name] = json.loads((root/'receipt.json').read_text())
        bindings[str(root/'receipt.json')] = sha(root/'receipt.json')
        for file, digest in receipts[name]['artifacts'].items():
            assert sha(root/file) == digest
            bindings[str(root/file)] = digest
    assert receipts['audit']['status'] == 'complete_refinement_disposition_and_diagnostic_readback'
    assert receipts['audit']['source_receipt_sha256'] == sha(roots['refinement']/'receipt.json')
    assert receipts['original']['status'] == 'complete_full_working_model_grid_export'
    original = pd.read_parquet(roots['original']/'unique_fits.parquet')
    full = pd.read_parquet(roots['original']/'full_settings.parquet')
    assert len(original) == 144040 and len(full) == 414720
    entries = [json.loads(line) for line in (roots['refinement']/'manifest.jsonl').read_text().splitlines()]
    assert len(entries) == receipts['audit']['dispositions'] == 658
    refinements = []
    for entry in entries:
        path = Path(entry['path'])
        assert path.resolve().is_relative_to(roots['refinement'].resolve()) and sha(path) == entry['sha256']
        bindings[str(path)] = entry['sha256']
        refinements.append((entry, json.loads(path.read_text())))
    updated, seen = overlay(original, refinements)
    expected = set(original.loc[original.numerical_review_required, KEY].itertuples(index=False, name=None))
    assert seen == expected and len(seen) == 658
    # Original columns are immutable; the new layer is additive.
    pd.testing.assert_frame_equal(updated[original.columns], original)
    extra = [column for column in updated if column not in original]
    expanded = full.merge(updated[KEY+extra], on=KEY, how='left', validate='many_to_one', sort=False)
    pd.testing.assert_frame_equal(expanded[full.columns], full)
    assert len(expanded) == 414720 and expanded.selection.notna().all()
    assert int(expanded.selection.ne('original_not_targeted').sum()) == 2099
    out = Path(config['output']); out.mkdir(parents=True, exist_ok=False)
    for name, frame in [('unique_fits.parquet', updated), ('full_settings.parquet', expanded)]:
        frame.to_parquet(out/name, index=False)
        pd.testing.assert_frame_equal(pd.read_parquet(out/name), frame)
    verify()
    for path, digest in bindings.items():
        assert sha(path) == digest
    write_json(out/'receipt.json', dict(status='complete_refined_estimate_overlay', plan_sha256=plan_hash,
        unique_fits=len(updated), full_setting_fits=len(expanded), refined_dispositions=len(seen),
        selected_review_unique_fits=int(updated.selected_review_required.sum()),
        selected_review_setting_fits=int(expanded.selected_review_required.sum()),
        selections=updated.selection.value_counts().to_dict(), source_bindings=bindings,
        artifacts={p.name: sha(p) for p in out.iterdir()},
        scope='Original columns preserved. Selected estimates use audited refinement candidates with all review flags retained; refinement errors retain originals. Conditional variances are not calibrated uncertainty; no biological effect or significance claim.'))


if __name__ == '__main__':
    main()
