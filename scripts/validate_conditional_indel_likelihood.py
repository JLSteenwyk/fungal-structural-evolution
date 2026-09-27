#!/usr/bin/env python3
"""Check ascertainment likelihood against independent matrix-exponential pruning."""
import json
from io import StringIO
from pathlib import Path
import numpy as np
from Bio import Phylo, SeqIO
from scipy.special import gammainc
from scipy.stats import gamma
from fit_conditional_indel_models import Likelihood
from replay_fastml_indel_probabilities import infer, analytic_check
from prepare_case_ancestral_neighborhoods import sha


def main():
    analytic_check()
    tree = Phylo.read(StringIO('(a:0.2,(b:0.3,c:0.4)I:0.1)R;'),'newick')
    seq = {'a':'01?1','b':'1?0?','c':'001?'}
    errors=[]
    for correction in ['all_taxa','observed_mask']:
        model=Likelihood(tree,seq,correction)
        for g,l,a in [(0.4,.9,.5),(2.,.02,2.),(.2,3.,.005),(.7,.5,100.)]:
            rates=4*np.diff(gammainc(a+1,a*gamma.ppf([0,.25,.5,.75,1],a=a,scale=1/a)))
            # Evaluate all original columns and each associated exclusion pattern explicitly.
            augmented={name:s+(''.join('?' if x=='?' else '0' for x in s) if correction=='observed_mask' else '0'*len(s)) for name,s in seq.items()}
            _,logp=infer(tree,augmented,g,l,rates)
            expected=float((logp[:4]-np.log(-np.expm1(logp[4:]))).sum())
            value=model.evaluate(np.log([g,l,a]))
            error=abs(expected-value)
            assert error<1e-10,(correction,g,l,a,error)
            errors.append(error)
    snapshot_path=Path('metadata/fastml_indel_probability_replay_snapshot_20260927.json')
    snapshot=json.loads(snapshot_path.read_text())
    comparisons=[]
    for previous in snapshot['jobs']:
        root=Path('results/ancestral/fastml-indels-20260927-v1')/previous['job_id']
        rp=root/'receipt.json'
        assert sha(rp)==previous['source_receipt_sha256']
        job=json.loads(rp.read_text())['job']
        tree=Phylo.read(root/'RESULTS/TheTree.INodes.ph','newick')
        seq={r.id:str(r.seq) for r in SeqIO.parse(job['characters'],'fasta')}
        value=Likelihood(tree,seq,'all_taxa').evaluate(np.log([previous['gain'],previous['loss'],previous['alpha']]))
        error=abs(value-previous['corrected_log_likelihood'])
        assert error<1e-8,(previous['job_id'],error)
        comparisons.append(dict(job_id=previous['job_id'],absolute_difference=error))
    result=dict(status='binary_working_likelihood_independent_checks_passed',
        analytic_enumeration='passed',synthetic_correction_comparisons=8,
        maximum_synthetic_error=max(errors),serialized_case_comparisons=comparisons,
        maximum_serialized_replay_difference=max(r['absolute_difference'] for r in comparisons),
        snapshot_sha256=sha(snapshot_path),script_sha256=sha(__file__),
        engine_sha256=sha('scripts/fit_conditional_indel_models.py'),
        scope='Numerical likelihood checks, including changed parameters, masks, gamma extremes and all81 existing independent case replays; not fitted-model validation, convergence or adequacy.')
    Path('metadata/conditional_indel_likelihood_validation_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],len(comparisons),result['maximum_serialized_replay_difference'])


if __name__=='__main__':
    main()
