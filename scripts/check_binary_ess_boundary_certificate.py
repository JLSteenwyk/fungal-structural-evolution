#!/usr/bin/env python3
"""Check exact zero-pair review certificates without relaxing numeric tolerance."""
import argparse
import json
from pathlib import Path

import arviz
import numpy as np

from ancestral_chain_diagnostics import diagnose as oracle
from independent_ancestral_scalar_diagnostics_v2 import diagnose, split
from binary_ess_boundary_certificate import candidates,compare_with_binary_certificate
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assert arviz.__version__=='0.22.0' and np.__version__=='2.2.6'
    p=Path('metadata/binary_ess_boundary_regression_20261002.json')
    regression=json.loads(p.read_text())
    pattern=np.asarray(regression['values'],dtype=np.uint8)
    assert pattern.shape==(4,50)
    values=(pattern==regression['state_index']).astype(float)
    original,separate=oracle(values),diagnose(values)
    errors,flags=compare_with_binary_certificate(original,separate,values)
    assert len(flags)==1 and flags[0]['metric']=='bulk_ess'
    exact=candidates(split(values));assert exact['exact_zero_pair_indices']==[1]
    assert any(abs(e-188.67018722995684)<1e-8 for e in exact['ess'])
    assert any(abs(e-199.59370238699853)<1e-8 for e in exact['ess'])
    rejected=[]
    for metric in ['rhat','bulk_ess','tail_ess','mean_mcse']:
        changed=dict(separate);changed[metric]+=1
        try:compare_with_binary_certificate(original,changed,values)
        except AssertionError:rejected.append(metric)
        else:raise AssertionError('Noncertified altered value accepted: '+metric)
    rng=np.random.default_rng(426);reviews=ordinary=0
    for draws in [20,21,50,75,500]:
        for frequency in [.01,.05,.5,.95,.99]:
            for repeat in range(10):
                values=(rng.random((4,draws))<frequency).astype(float)
                reference,separate=oracle(values),diagnose(values)
                errors,flags=compare_with_binary_certificate(reference,separate,values)
                reviews+=bool(flags);ordinary+=not bool(flags)
    paths=[__file__,'scripts/binary_ess_boundary_certificate.py','scripts/independent_ancestral_scalar_diagnostics_v2.py','scripts/ancestral_chain_diagnostics.py',str(p)]
    result=dict(status='passed_exact_binary_ess_zero_pair_review_contracts',
        actual_regression_certified_as_unresolved_not_agreement=True,real_case_certificate=exact,
        random_binary_fixtures=250,ordinary_agreement_fixtures=ordinary,numerical_review_fixtures=reviews,
        noncertified_metric_alterations_rejected=rejected,comparison_tolerance=dict(atol=1e-8,rtol=1e-8),
        source_hashes={str(p):sha(p) for p in paths},scientific_eligibility=False,
        production_integration_complete=False,scope='Exact binary integer/Fraction reachable-zero Geyer branch certificate. '
        'Differing values are explicitly unresolved, not passed numerical agreement; unrelated mismatches fail. '
        'No production stage restart or posterior/model qualification.')
    with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
