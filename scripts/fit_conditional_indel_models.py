#!/usr/bin/env python3
"""Fit binary gap working models with freshly evaluated ascertainment terms."""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from Bio import Phylo, SeqIO
from scipy.optimize import minimize
from scipy.special import gammainc, logsumexp
from scipy.stats import gamma
from threadpoolctl import threadpool_limits
from prepare_case_ancestral_neighborhoods import sha


class Likelihood:
    def __init__(self, tree, sequences, correction):
        assert correction in ['all_taxa', 'observed_mask']
        self.nodes = list(tree.find_clades(order='postorder'))
        self.index = {n:i for i,n in enumerate(self.nodes)}
        self.children = [[self.index[c] for c in n.clades] for n in self.nodes]
        self.lengths = [float(n.branch_length or 0) for n in self.nodes]
        assert min(self.lengths) >= 0
        tips = tree.get_terminals()
        assert set(sequences) == {n.name for n in tips}
        names = [n.name for n in tips]
        original = np.array([list(sequences[n]) for n in names]).T
        assert original.ndim == 2 and original.shape[0] > 0
        assert set(original.flatten()) <= {'0','1','?'}
        assert np.all((original == '1').any(axis=1))
        self.patterns, self.counts = np.unique(original, axis=0, return_counts=True)
        excluded = np.full(self.patterns.shape, '0')
        if correction == 'observed_mask':
            excluded[self.patterns == '?'] = '?'
        allpatterns, inverse = np.unique(np.concatenate([self.patterns,excluded]),axis=0,return_inverse=True)
        self.observed = inverse[:len(self.patterns)]
        self.excluded = inverse[len(self.patterns):]
        self.width = len(allpatterns)
        self.evidence = {}
        for i,n in enumerate(tips):
            letters = allpatterns[:,i]
            self.evidence[self.index[n]] = np.array([(letters=='?')|(letters==str(s)) for s in [0,1]],dtype=float)

    def evaluate(self, parameters):
        gain,loss,alpha = np.exp(parameters)
        pi = np.array([loss,gain])/(gain+loss)
        rates = 4*np.diff(gammainc(alpha+1, alpha*gamma.ppf([0,.25,.5,.75,1],a=alpha,scale=1/alpha)))
        assert np.all(rates > 0) and abs(rates.sum()-4) < 1e-10
        rate_logs = []
        for rate in rates:
            values,scales = {},{}
            for i,children in enumerate(self.children):
                part = self.evidence[i].copy() if i in self.evidence else np.ones((2,self.width))
                scale = np.zeros(self.width)
                for child in children:
                    # Stable analytic binary CTMC transition matrix, independent of expm replay.
                    changed = -np.expm1(-(gain+loss)*rate*self.lengths[child])
                    p = np.array([[1-pi[1]*changed,pi[1]*changed],
                                  [pi[0]*changed,1-pi[0]*changed]])
                    part *= p @ values.pop(child)
                    scale += scales.pop(child)
                    norm = part.max(axis=0)
                    if np.any(norm <= 0):
                        return -np.inf
                    part /= norm
                    scale += np.log(norm)
                values[i],scales[i] = part,scale
            root = len(self.nodes)-1
            rate_logs.append(np.log(pi @ values[root])+scales[root])
        logs = logsumexp(rate_logs,axis=0)-math.log(4)
        exclusion = logs[self.excluded]
        if np.any(exclusion >= 0):
            return -np.inf
        conditional = logs[self.observed]-np.log(-np.expm1(exclusion))
        return float(self.counts @ conditional)


def fit(job):
    with threadpool_limits(limits=1):
        return fit_serial(job)


def fit_serial(job):
    folder = Path(job['output'])
    folder.mkdir(parents=True,exist_ok=False)
    for p,h in job['pins'].items():
        assert sha(p)==h,p
    seq = {r.id:str(r.seq) for r in SeqIO.parse(job['characters'],'fasta')}
    result = dict(job=job)
    if not job['character_count']:
        result['status'] = 'no_coded_characters_no_inference'
    else:
        tree = Phylo.read(job['tree'],'newick')
        model = Likelihood(tree,seq,job['correction'])
        records = []
        bounds = [(-12,12),(-12,12),(math.log(.005),math.log(100))]
        for start in [(.2,.8,.5),(1.,1.,2.),(.02,2.,.1)]:
            x = np.log(start)
            def objective(v):
                value = model.evaluate(v)
                return -value if math.isfinite(value) else 1e100
            initial = model.evaluate(x)
            fitted = minimize(objective,x,method='L-BFGS-B',bounds=bounds,
                              options=dict(ftol=1e-11,gtol=1e-6,maxiter=500,maxls=40))
            value = model.evaluate(fitted.x)
            records.append(dict(start=list(start),initial_log_likelihood=initial,
                log_parameters=fitted.x.tolist(),parameters=np.exp(fitted.x).tolist(),
                log_likelihood=value,optimizer_success=bool(fitted.success),
                optimizer_message=str(fitted.message),iterations=int(fitted.nit),evaluations=int(fitted.nfev),
                gradient=fitted.jac.tolist(),
                near_bound=[bool(min(abs(v-lo),abs(v-hi))<1e-4) for v,(lo,hi) in zip(fitted.x,bounds)]))
        best = max(records,key=lambda r:r['log_likelihood'])
        result.update(status='three_start_fit_produced_pending_independent_validation',
            starts=records,best=best,log_parameter_bounds=bounds,
            fitted_log_likelihood_spread=max(r['log_likelihood'] for r in records)-min(r['log_likelihood'] for r in records))
    (folder/'receipt.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return dict(job_id=job['job_id'],status=result['status'],receipt_sha256=sha(folder/'receipt.json'))


def main():
    pp = Path('metadata/conditional_indel_model_plan_20260927.json')
    plan = json.loads(pp.read_text())
    for p,h in plan['pins'].items():
        assert sha(p)==h,p
    output = Path(plan['output'])
    output.mkdir(parents=True,exist_ok=False)
    with ProcessPoolExecutor(max_workers=2) as pool:
        results = []
        for result in pool.map(fit,plan['jobs']):
            results.append(result)
            print(result['job_id'],result['status'],flush=True)
    receipt = dict(status='all_working_model_fits_produced_pending_validation',
                   plan_sha256=sha(pp),jobs=results,scope=plan['scope'])
    (output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
