#!/usr/bin/env python3
"""Replay the diagnosed ML case through the guarded local-polish implementation."""
import json
from pathlib import Path
import numpy as np
from interior_gradient_polish import propose
from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_mixed_covariance import MatchedCovariance
from matched_ml_gradient import evaluate_gradient
from screen_duplication_alignment_reuse import sha


def main():
    pp = Path('metadata/refinement_gradient_diagnostic_plan_20260928.json')
    plan = json.loads(pp.read_text())
    bindings = {str(pp): sha(pp), **plan['pins']}
    for path, digest in bindings.items():
        assert sha(path) == digest
    original = json.loads(Path(plan['fit']).read_text())['payload']
    assert original['checks'] == dict(best_optimizer_success=True, projected_gradient_pass=False,
                                     upper_bound_contact=False, all_full_face_starts_agree=True,
                                     original_objective_not_worsened=True)
    with np.load(plan['inputs']) as a:
        values, bg, family, indices = a['matrix'], a['background'], a['family'], a['pattern_rows']
    with np.load(plan['factor']) as a:
        factor = a['factor'][indices]
    covariates = values[:, 1:][:, original['active_covariates']]
    x = np.column_stack([np.ones(len(values)), covariates / covariates.std(axis=0)])
    y = values[:, 0]
    cache = CachedMatchedML(bg, family, factor, x, y)
    def objective(t):
        return profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(t)), x, y)['negative_profiled_ml']
    def gradient(t):
        return evaluate_gradient(cache, np.expm1(t))['log1p_ratio_gradient']
    result = propose(objective, gradient, original['log1p_ratios'], np.log1p(original['maximum_ratio']))
    assert result['status'] == 'two_local_proposals_passed'
    independent = json.loads(Path('metadata/refinement_curvature_proposals_checked_20260928.json').read_text())
    for trial in result['proposals']:
        reference = next(c for c in independent['checks'] if c['hessian_step'] == trial['hessian_step'])
        np.testing.assert_allclose(trial['theta'], reference['candidate_theta'], rtol=0, atol=1e-14)
        np.testing.assert_allclose(trial['objective'], reference['objective'], rtol=0, atol=1e-9)
    for path, digest in bindings.items():
        assert sha(path) == digest
    proof = dict(status='passed_guarded_interior_polish_real_case_replay', source_hashes=bindings,
                 checker_sha256=sha(__file__), implementation_sha256=sha('scripts/interior_gradient_polish.py'),
                 separate_proposal_check_sha256=sha('metadata/refinement_curvature_proposals_checked_20260928.json'),
                 result=result, scope='One previously diagnosed case reproduced by guarded reusable routine; no stored fit or review flag changed. Wider applicability and integration require completed refinement audit.')
    with Path('metadata/interior_gradient_polish_case_checked_20260929.json').open('x') as handle:
        json.dump(proof, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(status=proof['status'], proposal_status=result['status'])))


if __name__ == '__main__':
    main()
