#!/usr/bin/env python3
"""Recompute fixed domain likelihoods by scaled pruning outside IQ-TREE."""
import csv,gzip,json,math,re
from io import StringIO
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,expm
from scipy.special import gammainc,logsumexp
from scipy.stats import gamma
from Bio import Phylo,SeqIO
from threadpoolctl import threadpool_limits
from prepare_case_ancestral_neighborhoods import sha

AA='ARNDCQEGHILKMFPSTWYV'


def matrices(path):
    text=Path(path).read_text();result={}
    for name in ['LG','WAG','JTT']:
        match=re.search(r'model '+name+r'=\s*(.*?);',text,re.S);assert match
        numbers=[float(x) for x in match.group(1).split()];assert len(numbers)==210
        matrix=np.zeros((20,20));k=0
        for i in range(1,20):
            for j in range(i):matrix[i,j]=matrix[j,i]=numbers[k];k+=1
        assert np.all(matrix[np.tril_indices(20,-1)]>0)
        result[name]=matrix
    return result


def generator(exchange,pi):
    q=exchange*pi[None,:];np.fill_diagonal(q,-q.sum(axis=1));q/=-(pi*np.diag(q)).sum()
    assert np.max(abs(q.sum(axis=1)))<1e-12 and np.max(abs(pi@q))<1e-12
    return q


def likelihood(tree,sequences,pi,q,rates):
    names=list(sequences);array=np.array([list(sequences[g]) for g in names]).T
    patterns,counts=np.unique(array,axis=0,return_counts=True)
    tipdata={g:np.array([[float(x=='-' or x==a) for x in patterns[:,i]] for a in AA]) for i,g in enumerate(names)}
    rootpi=np.sqrt(pi);symmetric=rootpi[:,None]*q/rootpi[None,:]
    assert np.max(abs(symmetric-symmetric.T))<1e-12
    values,vectors=eigh(symmetric);logs=[]
    for rate in rates:
        cache={};scales={}
        for node in tree.find_clades(order='postorder'):
            if not node.clades:
                cache[node]=tipdata[node.name];scales[node]=np.zeros(len(patterns));continue
            partial=np.ones((20,len(patterns)));scale=np.zeros(len(patterns))
            for child in node.clades:
                t=child.branch_length*rate
                if t<1e-8:p=expm(q*t)
                else:p=((vectors*np.exp(values*t))@vectors.T)*rootpi[None,:]/rootpi[:,None]
                assert p.min()>-1e-12 and np.max(abs(p.sum(axis=1)-1))<1e-10
                p=np.maximum(p,0.)
                partial*=p@cache.pop(child);scale+=scales.pop(child)
                norm=partial.max(axis=0);assert (norm>0).all();partial/=norm;scale+=np.log(norm)
            cache[node]=partial;scales[node]=scale
        logs.append(np.log(pi@cache[tree.root])+scales[tree.root])
    site=logsumexp(np.array(logs),axis=0)-math.log(len(rates))
    return float(site@counts)


def main():
    pp=Path('metadata/ancestral_domain_likelihood_replay_plan_20260927.json');plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():assert sha(p)==h
    source=Path(plan['fits']);rp=source/'receipt.json';r=json.loads(rp.read_text());assert sha(rp)==plan['fit_receipt_sha256']
    audit=json.loads(Path(plan['readback_receipt']).read_text());assert audit['source_receipt_sha256']==sha(rp)
    models=matrices(plan['exchangeability_source'])
    # Exact two-tip marginal likelihood identity tests both mixture weighting
    # and handling of an unobserved/gap state, using direct matrix exponentials.
    pi=np.arange(1,21,dtype=float);pi/=pi.sum();q=generator(models['LG'],pi)
    toy=Phylo.read(StringIO('(a:0.07,b:0.13);'),'newick');seq={'a':'AR-','b':'NAW'};rates=[.2,.6,1.1,2.1]
    expected=0.
    for a,b in zip(seq['a'],seq['b']):
        terms=[]
        for rate in rates:
            if a=='-':term=pi[AA.index(b)]
            else:term=pi[AA.index(a)]*expm(q*(.2*rate))[AA.index(a),AA.index(b)]
            terms.append(term)
        expected+=math.log(sum(terms)/4)
    assert abs(likelihood(toy,seq,pi,q,rates)-expected)<1e-11
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);rows=[]
    with (out/'likelihood_readback.tsv').open('w') as handle:
        writer=None
        for result in r['results']:
            job=result['job'];folder=source/job['job_id']
            assert sha(folder/'fit.ckp.gz')==result['artifacts']['fit.ckp.gz']
            assert sha(job['alignment'])==job['alignment_sha256']
            records=list(SeqIO.parse(job['alignment'],'fasta'));sequences={x.id:str(x.seq) for x in records};text=gzip.open(folder/'fit.ckp.gz','rt').read()
            tree=Phylo.read(StringIO(re.search(r'^ newick: (.+)$',text,re.M).group(1)),'newick')
            for tip in tree.get_terminals():tip.name=records[int(tip.name)].id
            assert {t.name for t in tree.get_terminals()}==set(sequences)
            alpha=float(re.search(r'^ gamma_shape: (.+)$',text,re.M).group(1));reference=float(re.search(r'^ 0: ([-+0-9.eE]+)',text,re.M).group(1))
            chars=''.join(sequences.values());assert set(chars)<=set(AA+'-')
            counts=np.array([chars.count(a) for a in AA]);pi=np.full(20,.05)
            for iteration in range(8):pi=(counts+chars.count('-')*pi)/len(chars)
            q=generator(models[job['model'].split('+')[0]],pi)
            bounds=gamma.ppf(np.linspace(0,1,5),a=alpha,scale=1/alpha);rates=4*np.diff(gammainc(alpha+1,alpha*bounds))
            value=likelihood(tree,sequences,pi,q,rates);error=value-reference
            row=dict(job_id=job['job_id'],checkpoint_log_likelihood=reference,replayed_log_likelihood=value,difference=error,passed=int(abs(error)<=plan['absolute_tolerance']))
            if writer is None:writer=csv.DictWriter(handle,list(row),delimiter='\t',lineterminator='\n');writer.writeheader()
            writer.writerow(row);handle.flush();rows.append(row);print(job['job_id'],error,flush=True)
    assert len(rows)==156
    result=dict(status='passed_all_156_independent_domain_likelihood_replays' if all(x['passed'] for x in rows) else 'domain_likelihood_replay_disagreements_require_review',fits=156,passing=sum(x['passed'] for x in rows),maximum_absolute_difference=max(abs(x['difference']) for x in rows),absolute_tolerance=plan['absolute_tolerance'],plan_sha256=sha(pp),source_receipt_sha256=sha(rp),analytic_two_tip_check='passed',script_sha256=sha(__file__),artifacts={'likelihood_readback.tsv':sha(out/'likelihood_readback.tsv')},scope='Independent scaled pruning and SciPy matrix exponentials/eigendecomposition; same empirical exchangeabilities, source alignment/tree and version-specific frequencies. Gaps marginalized as unknown, not reconstructed deletions. Reproducing saved likelihoods does not establish optimization convergence, model adequacy or ancestral states.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
