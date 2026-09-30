#!/usr/bin/env python3
"""Replay full-grid follow-up arithmetic against five verified saved cases."""
import argparse
import copy
import json
from pathlib import Path
from readback_whole_protein_flag_followup import check_result
from whole_protein_flag_followup import arrays
from screen_duplication_alignment_reuse import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text())
    bindings={str(args.plan):sha(args.plan),**plan['pins']}
    registry=json.loads(Path(plan['registry']).read_text())
    assert registry['status']=='completed_frozen_five_fit_numerical_recovery_candidates'
    assert registry['fits']==len(registry['candidates'])==5
    for path,h in registry['source_hashes'].items():
        assert path not in bindings or bindings[path]==h;bindings[path]=h
    def verify():
        for path,h in bindings.items():assert sha(path)==h,path
    verify();fit_plan=json.loads(Path(plan['fit_plan']).read_text())
    wanted={r['fit_input_id'] for r in registry['candidates']};recipes={}
    with (Path(fit_plan['inputs'])/'input_manifest.jsonl').open() as f:
        for line in f:
            row=json.loads(line)
            if row['fit_input_id'] in wanted:
                assert row['fit_input_id'] not in recipes;recipes[row['fit_input_id']]=row
    assert set(recipes)==wanted
    results=[];altered_rejected=False
    for entry in registry['candidates']:
        item=dict(fit_input_id=entry['fit_input_id'],tree=entry['tree'],
            path=entry['original_production_fit'],sha256=entry['original_production_sha256'])
        candidate=entry['candidate'];assert sha(candidate['path'])==candidate['sha256']
        saved=json.loads(Path(candidate['path']).read_text())
        if candidate['kind']=='audited_full_face_start_recovery':
            recovery=saved;refined=json.loads(Path(saved['source']).read_text())['payload']
            fitted=saved['fitted'];raw_beta=saved['raw_unit_beta'];raw_cov=saved['raw_unit_conditional_beta_covariance']
        else:
            assert candidate['kind']=='audited_all_face_refinement'
            recovery=None;refined=saved['payload'];fitted=refined
            raw_beta=refined['raw_unit_beta'];raw_cov=refined['raw_unit_conditional_beta_covariance']
        original,matrix,bg,family,factor,x,y,scales=arrays(item,recipes[item['fit_input_id']],fit_plan)
        result=dict(status='numerical_followup_passed_pending_independent_readback',refinement=refined,
            recovery=recovery,selected_theta=candidate['theta'],
            fitted={k:fitted[k] for k in ['negative_profiled_ml','beta','conditional_beta_covariance',
                'residual_quadratic','profiled_scale','variance_profile_denominator','residual_degrees_of_freedom']},
            raw_unit_beta=raw_beta,raw_unit_conditional_beta_covariance=raw_cov)
        checked=check_result(result,original,bg,family,factor,x,y,scales)
        assert checked['numerical_result_verified'] and checked['status']==result['status']
        results.append(dict(fit_input_id=item['fit_input_id'],tree=item['tree'],**checked))
        if not altered_rejected:
            altered=copy.deepcopy(result);altered['refinement']['candidates'][0]['objective']-=.1
            try:check_result(altered,original,bg,family,factor,x,y,scales)
            except AssertionError:altered_rejected=True
            else:raise AssertionError('Changed candidate likelihood accepted')
        print('Verified saved follow-up regression',len(results),'/5',flush=True)
    assert altered_rejected;verify()
    proof=dict(status='passed_five_saved_case_followup_arithmetic_regression',fits=5,
        likelihood_alteration_rejected=True,results=results,source_hashes=bindings,
        scope='All five previously verified saved numerical cases checked by the new full-grid checker. No optimizer rerun, no full-grid scientific result, and no substitute for the complete 375350-fit audit and all later flag follow-ups.')
    with Path(plan['output']).open('x') as f:json.dump(proof,f,indent=2);f.write('\n')


if __name__=='__main__':main()
