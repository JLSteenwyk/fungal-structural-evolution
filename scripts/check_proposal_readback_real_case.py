#!/usr/bin/env python3
"""Validate the separate proposal checker against the saved real ML case."""
import json
from pathlib import Path
import numpy as np
from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_mixed_covariance import MatchedCovariance
from matched_ml_gradient import evaluate_gradient
from readback_audited_gradient_proposals import check_proposal
from screen_duplication_alignment_reuse import sha


def main():
    source = Path('metadata/interior_gradient_polish_case_checked_20260929.json')
    saved = json.loads(source.read_text())
    plan = json.loads(Path('metadata/refinement_gradient_diagnostic_plan_20260928.json').read_text())
    bindings = {**saved['source_hashes'], str(source): sha(source),
                'scripts/readback_audited_gradient_proposals.py': sha('scripts/readback_audited_gradient_proposals.py')}
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    original = json.loads(Path(plan['fit']).read_text())['payload']
    with np.load(plan['inputs']) as arrays:
        values, bg, family, rows = arrays['matrix'], arrays['background'], arrays['family'], arrays['pattern_rows']
    with np.load(plan['factor']) as arrays:
        factor = arrays['factor'][rows]
    active = np.ptp(values[:, 1:], axis=0) > 1e-12
    scales = np.std(values[:, 1:][:, active], axis=0)
    np.testing.assert_array_equal(active, original['active_covariates'])
    np.testing.assert_array_equal(scales, original['covariate_scales'])
    x = np.column_stack((np.ones(len(values)), values[:, 1:][:, active] / scales))
    y = values[:, 0]
    cache = CachedMatchedML(bg, family, factor, x, y)
    evaluations = dict(objective=0, gradient=0)
    def objective(theta):
        evaluations['objective'] += 1
        return profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(theta)), x, y)['negative_profiled_ml']
    def gradient(theta):
        evaluations['gradient'] += 1
        return evaluate_gradient(cache, np.expm1(theta))['log1p_ratio_gradient']
    check_proposal(saved['result'], objective, gradient, original['log1p_ratios'], np.log1p(original['maximum_ratio']))
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    proof = dict(status='passed_separate_proposal_readback_real_case', records=len(y), proposals=2,
                 evaluations=evaluations, source_hashes=bindings, script_sha256=sha(__file__),
                 scope='Saved real-case proposals replayed without correction helper; shared likelihood and gradient libraries. No replacement or full-grid acceptance.')
    with Path('metadata/proposal_readback_real_case_completed_20260929.json').open('x') as h:
        json.dump(proof, h, indent=2)
        h.write('\n')
    print(json.dumps({k:v for k,v in proof.items() if k != 'source_hashes'}), flush=True)


if __name__ == '__main__':
    main()
