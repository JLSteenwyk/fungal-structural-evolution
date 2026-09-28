"""Exercise source joins and unchanged-grid invariants on completed refinements."""
import copy
import json
from pathlib import Path

import pandas as pd
from ancestral_chain_attempt import sha, write_json
from integrate_matched_refinements import overlay, KEY


def main():
    original = pd.read_parquet('results/structural_comparisons/full-working-model-grid-export-20260927-v1/unique_fits.parquet')
    root = Path('results/structural_comparisons/matched-reml-analytic-refinement-20260928-v1')
    out = Path('results/model_validation/matched-refinement-overlay-checks-20260928-v2')
    out.mkdir(parents=True, exist_ok=False)
    fixtures = []
    for path in sorted(root.glob('*.json')):
        saved = json.loads(path.read_text())
        if 'payload' not in saved:
            continue
        snapshot = out/path.name
        snapshot.write_bytes(path.read_bytes())
        entry = dict(fit_input_id=saved['source']['fit_input_id'], tree=saved['source']['tree'],
                     sha256=sha(snapshot), status=saved['payload']['status'])
        fixtures.append((entry, saved))
    assert fixtures
    result, seen = overlay(original, fixtures)
    pd.testing.assert_frame_equal(result[original.columns], original)
    assert len(seen) == len(fixtures)
    assert result.selection.ne('original_not_targeted').sum() == len(fixtures)
    full = pd.read_parquet('results/structural_comparisons/full-working-model-grid-export-20260927-v1/full_settings.parquet')
    extra = [c for c in result if c not in original]
    expanded = full.merge(result[KEY+extra], on=KEY, how='left', validate='many_to_one', sort=False)
    pd.testing.assert_frame_equal(expanded[full.columns], full)
    assert len(expanded) == 414720 and expanded.selection.notna().all()
    lookup = result.set_index(KEY)
    for entry, saved in fixtures:
        p = saved['payload']; row = lookup.loc[(entry['fit_input_id'], entry['tree'])]
        assert row.selected_intercept == p['raw_unit_beta'][0]
        assert row.selected_conditional_intercept_variance == p['raw_unit_conditional_beta_covariance'][0][0]
        assert row.selected_negative_profiled_reml == p['negative_profiled_reml']
    entry, saved = fixtures[0]
    small = original[(original.fit_input_id == entry['fit_input_id']) & (original.tree == entry['tree'])].reset_index(drop=True)
    failed = copy.deepcopy(saved); failed['payload'] = dict(status='refinement_error_requires_review')
    error_entry = {**entry, 'status': failed['payload']['status']}
    retained, _ = overlay(small, [(error_entry, failed)])
    assert retained.iloc[0].selected_intercept == small.iloc[0].intercept
    assert retained.iloc[0].selected_review_required
    assert retained.iloc[0].selection == 'original_after_refinement_error'
    bad_hash = copy.deepcopy(saved); bad_hash['source']['source_fit_sha256'] = 'changed'
    worse = copy.deepcopy(saved); worse['payload']['negative_profiled_reml'] = small.iloc[0].negative_profiled_reml + 1
    rejected = 0
    for records in [[(entry, saved), (entry, saved)], [(entry, bad_hash)], [(entry, worse)]]:
        try:
            overlay(small, records)
        except AssertionError:
            rejected += 1
        else:
            raise AssertionError('Invalid refinement accepted')
    write_json(out/'receipt.json', dict(status='passed_frozen_refinement_overlay_checks',
        frozen_completed_fixtures=len(fixtures), original_rows_preserved=len(original),
        invalid_cases_rejected=rejected, error_fallback_checked=True,
        artifacts={p.name: sha(p) for p in out.iterdir()},
        scope='Completed snapshot used to test export mapping; not full production completion or independent likelihood verification.'))
    write_json(Path('metadata/matched_refinement_overlay_checks_20260928.json'), dict(
        status='passed_frozen_refinement_overlay_checks', fixtures=len(fixtures),
        receipt=str(out/'receipt.json'), receipt_sha256=sha(out/'receipt.json'),
        script_sha256=sha(__file__), integration_script_sha256=sha('scripts/integrate_matched_refinements.py')))
    print('Passed', len(fixtures), 'frozen refinements and', rejected, 'rejection cases')


if __name__ == '__main__':
    main()
