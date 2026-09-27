#!/usr/bin/env python3
"""Validate direct ascertainment probabilities against independent replays."""
import json
from pathlib import Path
from io import StringIO
import numpy as np
from Bio import Phylo,SeqIO
from stable_indel_likelihood import StableLikelihood
from diagnose_indel_extreme_rates import high_precision
from prepare_case_ancestral_neighborhoods import sha


def main():
    tree=Phylo.read(StringIO('(a:0.2,(b:0.3,c:0.4)I:0.1)R;'),'newick')
    seq={'a':'01?1','b':'1?0?','c':'001?'}
    synthetic=[]
    for correction in ['all_taxa','observed_mask']:
        model=StableLikelihood(tree,seq,correction)
        for p in [[.4,.9,.5],[.00001,160000.,100.],[.2,3.,.005]]:
            expected=float(high_precision(tree,seq,p,correction))
            observed=model.evaluate(np.log(p))
            error=abs(observed-expected)
            assert error<1e-10,(correction,p,error)
            synthetic.append(dict(correction=correction,parameters=p,error=error))
    comparisons=[]
    root=Path('results/ancestral/conditional-indel-fit-audit-20260927-v2')
    for ap in sorted(root.glob('*.json')):
        if ap.name=='receipt.json':
            continue
        audit=json.loads(ap.read_text())
        if 'starts' not in audit:
            continue
        rp=Path('results/ancestral/conditional-indel-models-20260927-v1')/audit['job_id']/'receipt.json'
        assert sha(rp)==audit['source_receipt_sha256']
        fit=json.loads(rp.read_text());job=fit['job']
        seq={r.id:str(r.seq) for r in SeqIO.parse(job['characters'],'fasta')}
        model=StableLikelihood(Phylo.read(job['tree'],'newick'),seq,job['correction'])
        for start,reference in zip(fit['starts'],audit['starts']):
            observed=model.evaluate(start['log_parameters'])
            error=abs(observed-reference['independent_log_likelihood'])
            assert error<1e-6,(audit['job_id'],reference['start'],error)
            comparisons.append(dict(job_id=audit['job_id'],start=reference['start'],
                independent_precision=reference['precision'],error=error,
                source_audit_sha256=sha(ap)))
    result=dict(status='stable_event_likelihood_validation_snapshot_passed',
        synthetic_checks=synthetic,fit_comparisons=comparisons,
        maximum_error=max(r['error'] for r in comparisons),
        maximum_high_precision_error=max(r['error'] for r in comparisons if r['independent_precision']=='70_digit_recursive_pruning'),
        pins={p:sha(p) for p in ['scripts/stable_indel_likelihood.py','scripts/fit_conditional_indel_models.py','scripts/diagnose_indel_extreme_rates.py',__file__]},
        scope='Direct event-denominator arithmetic validated against independent calculations at existing parameter values. Snapshot includes all completed audit records discovered on entry. Not new optimization or completed ancestral inference.')
    Path('metadata/stable_indel_likelihood_validation_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    print(len(comparisons),result['maximum_error'],result['maximum_high_precision_error'])


if __name__=='__main__':
    main()
