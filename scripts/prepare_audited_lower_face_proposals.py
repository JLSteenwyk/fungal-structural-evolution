#!/usr/bin/env python3
"""Prepare guarded proposals only after complete refinement-output replay."""
import json
import subprocess
import time
from collections import Counter
from pathlib import Path
import numpy as np
import psutil
from lower_face_gradient_polish import propose
from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_mixed_covariance import MatchedCovariance
from matched_ml_gradient import evaluate_gradient
from screen_duplication_alignment_reuse import sha


def main():
    pp = Path('metadata/audited_lower_face_proposals_plan_20260929.json')
    plan = json.loads(pp.read_text())
    bindings = {str(pp): sha(pp), **plan['pins']}
    launch = json.loads(Path(plan['audit_launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            p = psutil.Process(launch['pid'])
            if abs(p.create_time() - launch['created']) > .01 or p.status() == psutil.STATUS_ZOMBIE:
                break
            assert p.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(x.split('=', 1) for x in subprocess.check_output(['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    root = Path(plan['refinements'])
    receipt = json.loads((root / 'receipt.json').read_text())
    proof = json.loads(Path(plan['audit_proof']).read_text())
    assert proof['status'] == 'passed_full_flagged_polynomial_refinement_replay'
    assert proof['source_receipt_sha256'] == sha(root / 'receipt.json') and proof['fits'] == 601
    bindings[plan['audit_proof']] = sha(plan['audit_proof'])
    bindings[str(root / 'receipt.json')] = sha(root / 'receipt.json')
    bindings.update(receipt['source_hashes'])
    bindings.update({str(root / n): h for n, h in receipt['artifacts'].items()})
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    inputs = Path(plan['inputs'])
    factors = {}
    for path in Path(plan['factors']).glob('*.npz'):
        with np.load(path) as arrays:
            factors[path.stem] = arrays['factor']
    manifest = json.loads((root / 'fit_manifest.json').read_text())
    assert len(manifest) == 601
    out = Path(plan['output'])
    out.mkdir(exist_ok=False)
    counts, seen = Counter(), set()
    eligible_checks = dict(best_optimizer_success=True, projected_gradient_pass=False,
                           upper_bound_contact=False, all_full_face_starts_agree=True,
                           original_objective_not_worsened=True)
    with (out / 'dispositions.jsonl').open('w') as handle:
        for row in manifest:
            key = row['fit_input_id'], row['tree']
            assert key not in seen
            seen.add(key)
            path = root / row['path']
            assert sha(path) == row['sha256']
            saved = json.loads(path.read_text())
            p = saved['payload']
            result = dict(fit_input_id=key[0], tree=key[1], source_fit=str(path), source_sha256=row['sha256'], original_status=p['status'])
            if p['status'] == 'ml_refinement_passed_numerical_checks':
                result['disposition'] = 'already_passed_no_proposal'
            elif p['checks'] != eligible_checks:
                result['disposition'] = 'other_review_flags_preserved'
                result['original_checks'] = p['checks']
            else:
                input_path = inputs / (key[0] + '.npz')
                assert sha(input_path) == saved['input_sha256']
                with np.load(input_path) as arrays:
                    values, bg, family = arrays['matrix'], arrays['background'], arrays['family']
                    factor = factors[key[1]][arrays['pattern_rows']]
                x = np.column_stack([np.ones(len(values)), values[:, 1:][:, p['active_covariates']] / p['covariate_scales']])
                y = values[:, 0]
                cache = CachedMatchedML(bg, family, factor, x, y)
                def objective(theta):
                    return profiled_ml(MatchedCovariance(bg, family, factor, 1., *np.expm1(theta)), x, y)['negative_profiled_ml']
                def gradient(theta):
                    return evaluate_gradient(cache, np.expm1(theta))['log1p_ratio_gradient']
                np.testing.assert_allclose(objective(np.asarray(p['log1p_ratios'])), p['negative_profiled_ml'], rtol=1e-10, atol=1e-8)
                result['proposal'] = propose(objective, gradient, p['log1p_ratios'], np.log1p(p['maximum_ratio']))
                result['disposition'] = result['proposal']['status']
            counts[result['disposition']] += 1
            handle.write(json.dumps(result, allow_nan=False, separators=(',', ':')) + '\n')
    verify()
    receipt = dict(status='complete_audited_lower_face_proposals_pending_independent_readback', fits=len(seen),
                   counts=dict(counts), source_hashes=bindings, script_sha256=sha(__file__),
                   artifacts={'dispositions.jsonl': sha(out / 'dispositions.jsonl')},
                   scope='All601 audited frozen refinements retained. Only gradient-only mixed lower-face cases receive local proposals; other cases retain explicit inapplicable dispositions. Existing fit statuses and all other review flags unchanged. No proposed point selected or integrated; full proposal readback and full-grid audit remain required.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'source_hashes'}), flush=True)


if __name__ == '__main__':
    main()
