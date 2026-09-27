#!/usr/bin/env python3
"""Evaluate finite likelihood neighborhoods, without claiming full optimization."""
import csv,gzip,json,re,itertools
from pathlib import Path
from io import StringIO
import numpy as np
from scipy.special import gammainc
from scipy.stats import gamma
from Bio import Phylo,SeqIO
from threadpoolctl import threadpool_limits
from prepare_case_ancestral_neighborhoods import sha
from replay_ancestral_domain_likelihoods import AA,matrices,generator,likelihood


def main():
    pp=Path('metadata/ancestral_domain_fit_neighborhood_plan_20260927.json');plan=json.loads(pp.read_text())
    def verify():
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();root=Path(plan['fits']);fits=json.loads((root/'receipt.json').read_text())['results'];params={x['job_id']:x for x in json.loads(Path(plan['parameters']).read_text())};models=matrices(plan['exchangeabilities'])
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);rows=[];summaries=[]
    for result in fits:
        job=result['job'];name=job['job_id'];ck=root/name/'fit.ckp.gz';assert sha(ck)==result['artifacts']['fit.ckp.gz'] and sha(job['alignment'])==job['alignment_sha256']
        records=list(SeqIO.parse(job['alignment'],'fasta'));seq={r.id:str(r.seq) for r in records};text=gzip.open(ck,'rt').read();tree=Phylo.read(StringIO(re.search(r'^ newick: (.+)$',text,re.M).group(1)),'newick')
        for tip in tree.get_terminals():tip.name=records[int(tip.name)].id
        original={n:n.branch_length for n in tree.find_clades() if n is not tree.root};par=params[name];alpha=par['gamma_shape'];pi=np.array([par['reconstructed_iqtree_empirical_frequencies'][a] for a in AA]);q=generator(models[job['model'].split('+')[0]],pi)
        expected=float(re.search(r'^ 0: ([-+0-9.eE]+)',text,re.M).group(1));base=likelihood(tree,seq,pi,q,par['gamma_category_rates']);assert abs(base-expected)<.001
        local=[]
        for af,bf in itertools.product(plan['alpha_factors'],plan['branch_factors']):
            a=alpha*af;bounds=gamma.ppf(np.linspace(0,1,5),a=a,scale=1/a);rates=4*np.diff(gammainc(a+1,a*bounds));assert np.isfinite(rates).all() and (rates>=0).all() and abs(rates.mean()-1)<1e-10
            for n,length in original.items():n.branch_length=length*bf
            value=likelihood(tree,seq,pi,q,rates);assert np.isfinite(value)
            if af==bf==1:assert abs(value-base)<1e-8
            row=dict(job_id=name,family=job['family'],boundary=job['boundary'],method=job['method'],model=job['model'],alpha_factor=af,branch_factor=bf,gamma_shape=a,below_original_alpha_lower_bound=int(a<plan['original_alpha_lower_bound']),log_likelihood=value,baseline_log_likelihood=base,improvement=value-base)
            rows.append(row);local.append(row)
        best=max(local,key=lambda x:x['improvement']);summaries.append(dict(job_id=name,baseline_gamma_shape=alpha,maximum_grid_improvement=best['improvement'],best_alpha_factor=best['alpha_factor'],best_branch_factor=best['branch_factor'],improved_above_tolerance=int(best['improvement']>plan['improvement_tolerance']),best_below_original_alpha_bound=best['below_original_alpha_lower_bound']))
        # Per-fit checkpoints are retained even if a subsequent fit fails.
        (out/(name+'.json')).write_text(json.dumps(dict(grid=local,summary=summaries[-1]),indent=2)+'\n');print(name,best['improvement'],best['alpha_factor'],best['branch_factor'],flush=True)
    assert len(summaries)==156 and len(rows)==156*len(plan['alpha_factors'])*len(plan['branch_factors'])
    for filename,data in [('grid.tsv',rows),('fit_summary.tsv',summaries)]:
        with (out/filename).open('w') as h:
            w=csv.DictWriter(h,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    verify();r=dict(status='complete_finite_domain_likelihood_neighborhoods',fits=156,grid_points=len(rows),fits_improved_above_tolerance=sum(x['improved_above_tolerance'] for x in summaries),maximum_grid_improvement=max(x['maximum_grid_improvement'] for x in summaries),plan_sha256=sha(pp),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Finite two-parameter diagnostic with fixed frequencies/topology and uniform branch scaling. An improving point disproves optimality in these directions; absence of improvement does not prove convergence, global optimality, model adequacy, or stable ancestral probabilities. Below-bound evaluations extend the parameter domain and do not alone imply failure within original constraints.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='artifacts'}),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=1):main()
