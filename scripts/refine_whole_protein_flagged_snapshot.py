#!/usr/bin/env python3
"""Refine an explicitly frozen, previously replayed set; never replace production fits."""
import argparse
import json
from pathlib import Path
import numpy as np
from refine_matched_ml_analytic import refine
from readback_selected_refinement_candidates import independent_fit
from matched_mixed_covariance import MatchedCovariance
from screen_duplication_alignment_reuse import sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    def verify():
        for p, h in bindings.items():
            assert sha(p) == h, p
    verify()
    proof = json.loads(Path(plan['original_replay']).read_text())
    assert proof['status'] == 'completed_frozen_prefix_whole_protein_flag_replay'
    bindings.update(proof['source_hashes'])
    verify()
    wanted = {(r['fit_input_id'], r['tree']) for r in proof['results']}
    rows = [json.loads(line) for line in Path(plan['frozen_manifest']).read_text().splitlines()]
    selected = [r for r in rows if (r['fit_input_id'], r['tree']) in wanted]
    assert len(selected) == len(wanted) == proof['flagged_dispositions']
    root = Path(plan['inputs'])
    input_receipt = json.loads((root/'receipt.json').read_text())
    assert sha(root/'input_manifest.jsonl') == input_receipt['artifacts']['input_manifest.jsonl']
    ids = {k[0] for k in wanted}
    inputs = {}
    for line in (root/'input_manifest.jsonl').open():
        r = json.loads(line)
        if r['fit_input_id'] in ids:
            assert r['fit_input_id'] not in inputs
            inputs[r['fit_input_id']] = r
    assert set(inputs) == ids
    output = Path(plan['output'])
    output.mkdir(parents=True, exist_ok=False)
    dispositions = []
    for item in selected:
        source = Path(item['path'])
        assert sha(source) == item['sha256'] == bindings[str(source)]
        original = json.loads(source.read_text())
        assert original['payload']['status'] == 'ml_candidate_requires_optimization_review'
        recipe = inputs[item['fit_input_id']]
        path = root/recipe['path']
        assert sha(path) == recipe['sha256'] == original['input_sha256']
        assert original['specification'] == recipe['recipe']['specification']
        with np.load(path, allow_pickle=False) as a:
            matrix, bg, family, index = a['matrix'], a['background'], a['family'], a['pattern_rows']
        factor_path = Path(plan['factors'])/(item['tree']+'.npz')
        assert sha(factor_path) == bindings[str(factor_path)]
        with np.load(factor_path, allow_pickle=False) as a:
            factor = a['factor'][index]
        scales = np.std(matrix[:, 2:], axis=0)
        np.testing.assert_array_equal(scales, original['payload']['covariate_scales'])
        assert np.all(matrix[:, 1] == 1)
        design = np.column_stack([matrix[:, 1], matrix[:, 2:]/scales])
        response = matrix[:, 0]
        result = refine(bg, family, factor, design, response,
                        original['payload']['log1p_ratios'], maximum_ratio=10000., maxiter=1000)
        assert len(result['candidates']) == 24
        np.testing.assert_allclose(result['original_objective'], original['payload']['negative_profiled_ml'], rtol=1e-9, atol=1e-7)
        max_error = 0.
        for candidate in result['candidates']:
            fit = independent_fit(MatchedCovariance(bg, family, factor, 1., *np.expm1(candidate['theta'])), design, response)
            error = abs(fit['negative_profiled_ml']-candidate['objective'])
            np.testing.assert_allclose(fit['negative_profiled_ml'], candidate['objective'], rtol=1e-9, atol=1e-7)
            max_error = max(max_error, error)
        best = min(result['candidates'], key=lambda c: (c['objective'], not c['success']))
        np.testing.assert_array_equal(best['theta'], result['log1p_ratios'])
        fit = independent_fit(MatchedCovariance(bg, family, factor, 1., *np.expm1(best['theta'])), design, response)
        for key, value in fit.items():
            np.testing.assert_allclose(result[key], value, rtol=1e-7, atol=1e-8)
        conversion = np.r_[1., 1/scales]
        result['raw_unit_beta'] = (fit['beta']*conversion).tolist()
        result['raw_unit_conditional_beta_covariance'] = (fit['conditional_beta_covariance']*np.outer(conversion, conversion)).tolist()
        record = dict(fit_input_id=item['fit_input_id'], tree=item['tree'], source_fit=str(source),
                      source_fit_sha256=item['sha256'], input_sha256=recipe['sha256'],
                      columns=original['specification']['columns'], covariate_scales=scales.tolist(),
                      payload=result, candidate_likelihoods_replayed=24,
                      maximum_independent_objective_error=max_error)
        target = output/(item['fit_input_id']+'-'+item['tree']+'.json')
        target.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
        dispositions.append(dict(fit_input_id=item['fit_input_id'], tree=item['tree'],
                                 path=target.name, sha256=sha(target), status=result['status'],
                                 checks=result['checks'], objective_improvement=result['objective_improvement']))
        print(json.dumps(dispositions[-1]), flush=True)
    verify()
    result = dict(status='complete_frozen_whole_protein_refinement_pending_serialized_review',
                  fits=len(dispositions), candidate_likelihoods_replayed=24*len(dispositions),
                  dispositions=dispositions, source_hashes=bindings,
                  scope='All flagged fits in the previously audited frozen prefix refined separately with unchanged numerical acceptance thresholds. All candidates replayed using independent GLS normal equations with the shared covariance operator. Original production outputs and flags unchanged. Serialized review, integration, later flags and calibrated inference remain outstanding.')
    (output/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
