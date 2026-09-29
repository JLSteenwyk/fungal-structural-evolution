#!/usr/bin/env python3
"""Exercise complete payload replay and deliberate corruption detection."""
import copy
import json
from pathlib import Path
import numpy as np
from fit_whole_protein_ml import fit_input
from check_whole_protein_fit_payload import check_payload
from screen_duplication_alignment_reuse import sha


def main():
    rng=np.random.default_rng(91823)
    n=36
    bg=np.repeat(np.arange(9),4)
    family=bg//3
    matrix=np.column_stack([rng.normal(size=n),np.ones(n),rng.normal(size=(n,4))])
    columns=['rmsd','intercept','a','b','c','d']
    results=[]
    for factor in [np.zeros((n,0)),rng.normal(size=(n,3))*.2]:
        payload=fit_input(bg,family,factor,matrix,columns)
        results.append(check_payload(payload,bg,family,factor,matrix,columns))
    rejected=[]
    for kind in ['objective','coefficient','covariance','gradient','status','candidate_omission']:
        bad=copy.deepcopy(payload)
        if kind=='objective':bad['candidates'][0]['objective']+=1
        elif kind=='coefficient':bad['raw_unit_beta'][0]+=1
        elif kind=='covariance':bad['raw_unit_conditional_beta_covariance'][0][0]+=1
        elif kind=='gradient':bad['projected_gradient'][0]+=1
        elif kind=='status':bad['checks']['projected_gradient_pass']=not bad['checks']['projected_gradient_pass']
        else:bad['candidates'].pop()
        try:check_payload(bad,bg,family,factor,matrix,columns)
        except AssertionError:rejected.append(kind)
        else:raise AssertionError(kind+' corruption accepted')
    result=dict(status='passed_whole_protein_payload_replay_fixtures',valid_cases=results,rejected_corruptions=rejected,source_hashes={p:sha(p) for p in ['scripts/check_whole_protein_fit_payload.py','scripts/check_whole_protein_payload_fixtures.py','scripts/fit_whole_protein_ml.py','scripts/readback_selected_refinement_candidates.py']},scope='Synthetic payload checks. Separate normal-equation likelihood implementation shares covariance solver and analytic-gradient library. No production audit, global-optimum proof, calibrated uncertainty or biological inference.')
    with Path('metadata/whole_protein_payload_replay_fixtures_20260929.json').open('x') as handle:
        json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
