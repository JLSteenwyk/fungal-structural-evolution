#!/usr/bin/env python3
"""Audit every marginal with fixed-root inside/outside messages and direct expm."""
import csv,gzip,json,re,hashlib,subprocess,time
import psutil
from map_ancestral_fitted_nodes import signatures
from io import StringIO
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.special import logsumexp
from Bio import Phylo,SeqIO
from threadpoolctl import threadpool_limits
from prepare_case_ancestral_neighborhoods import sha,read
from replay_ancestral_domain_likelihoods import AA,matrices


def independent(tree,sequences,pi,exchange,rates,targets):
    # Construct Q separately. No spectral decomposition, pattern compression,
    # traversal rerooting, or producer posterior routine is used.
    q=np.zeros((20,20))
    for i in range(20):
        for j in range(20):
            if i!=j:q[i,j]=exchange[i,j]*pi[j]
        q[i,i]=-q[i].sum()
    q/=sum(pi[i]*(-q[i,i]) for i in range(20))
    width=len(next(iter(sequences.values())))
    evidence={}
    for name,seq in sequences.items():
        a=np.ones((20,width))
        for k,letter in enumerate(seq):
            if letter!='-':a[:,k]=0;a[AA.index(letter),k]=1
        evidence[name]=a
    probabilities=[];loglikes=[]
    for rate in rates:
        transition={c:expm(q*(c.branch_length*rate)) for n in tree.find_clades() for c in n.clades}
        assert all(p.min()>=0 and np.max(abs(p.sum(axis=1)-1))<1e-10 for p in transition.values())
        inside={};scale={};message={}
        for n in tree.find_clades(order='postorder'):
            v=evidence[n.name].copy() if n.is_terminal() else np.ones((20,width))
            s=np.zeros(width)
            for c in n.clades:
                message[c]=transition[c]@inside[c];v*=message[c];s+=scale[c]
                z=v.sum(axis=0);assert (z>0).all();v/=z;s+=np.log(z)
            inside[n]=v;scale[n]=s
        loglikes.append(np.log(pi@inside[tree.root])+scale[tree.root])
        outside={tree.root:np.broadcast_to(pi[:,None],(20,width)).copy()}
        for n in tree.find_clades(order='preorder'):
            for c in n.clades:
                v=outside[n].copy()
                for sibling in n.clades:
                    if sibling is c:continue
                    v*=message[sibling];z=v.sum(axis=0);assert (z>0).all();v/=z
                v=transition[c].T@v;z=v.sum(axis=0);assert (z>0).all();outside[c]=v/z
        result=[]
        for n in targets:
            v=inside[n]*outside[n];v/=v.sum(axis=0);result.append(v.T)
        probabilities.append(result)
    ll=np.array(loglikes);weights=np.exp(ll-logsumexp(ll,axis=0))
    result=np.einsum('rnsa,rs->nsa',np.array(probabilities),weights)
    return result,logsumexp(ll,axis=0)-np.log(len(rates))


def analytic(exchange):
    pi=np.arange(1,21,dtype=float);pi/=pi.sum();rates=[.2,.6,1.1,2.1]
    t=Phylo.read(StringIO('((a:0.07,b:0.13)n1:0.09,c:0.11)n0;'),'newick')
    seq={'a':'A-','b':'NW','c':'R-'};targets=[t.root,t.root.clades[0]]
    got,ll=independent(t,seq,pi,exchange,rates,targets)
    q=exchange*pi[None,:];np.fill_diagonal(q,-q.sum(axis=1));q/=-(pi*np.diag(q)).sum()
    expected=np.zeros((2,2,20));likes=[]
    for col in range(2):
        joint=np.zeros((20,20))
        for rate in rates:
            edge=expm(q*.09*rate)
            a=expm(q*.07*rate)[:,AA.index('A')] if col==0 else np.ones(20)
            b=expm(q*.13*rate)[:,AA.index(seq['b'][col])]
            c=expm(q*.11*rate)[:,AA.index('R')] if col==0 else np.ones(20)
            joint+=(pi*c)[:,None]*edge*(a*b)[None,:]/4
        likes.append(np.log(joint.sum()));expected[0,col]=joint.sum(axis=1)/joint.sum();expected[1,col]=joint.sum(axis=0)/joint.sum()
    assert np.max(abs(expected-got))<1e-11 and np.max(abs(ll-likes))<1e-11


def main():
    pp=Path('metadata/whole_protein_posterior_audit_plan_20260927.json');plan=json.loads(pp.read_text())
    def verify():
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();identity=plan['producer']
    while True:
        state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',identity['unit'],'--property=ActiveState,MainPID,Result,ExecMainStatus'],text=True).splitlines())
        if state['ActiveState']=='inactive':
            assert state['Result']=='success' and state['ExecMainStatus']=='0';break
        assert state['ActiveState'] in ['active','activating','deactivating'],state
        proc=psutil.Process(identity['pid']);assert proc.create_time()==identity['created'] and proc.cmdline()==identity['cmdline'] and int(state['MainPID'])==proc.pid
        print('waiting_for_verified_whole_protein_marginals',flush=True);time.sleep(30)
    verify();source=Path(plan['source']);sr=json.loads((source/'receipt.json').read_text())
    assert sr['status']=='all_78_whole_protein_marginals_produced_pending_independent_audit' and sr['fits']==78 and sr['nodes']==234
    assert sr['plan_sha256']==sha(plan['producer_plan'])
    base_plan=json.loads(Path(plan['producer_plan']).read_text());ar=Path(base_plan['audit_receipt']);audit=json.loads(ar.read_text())
    assert sr['audit_receipt_sha256']==sha(ar)
    assert audit['source_receipt_sha256']==sha(Path(base_plan['fits'])/'receipt.json')
    assert audit['artifacts']['parameters.json']==sha(base_plan['parameters'])
    for p,h in sr['artifacts'].items():assert sha(source/p)==h
    base=json.loads(Path(plan['producer_plan']).read_text());models=matrices(base['exchangeabilities']);analytic(models['LG'])
    fits=Path(base['fits']);jobs=json.loads((fits/'receipt.json').read_text())['results'];params={p['job_id']:p for p in json.loads(Path(base['parameters']).read_text())}
    nodes=[r for r in read(source/'candidate_mapping.tsv') if r['status']=='matched_unrooted_vertex']
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);rows=[]
    for fit in jobs:
        job=fit['job'];name=job['job_id'];ck=fits/name/'fit.ckp.gz';assert sha(ck)==fit['artifacts']['fit.ckp.gz'];assert sha(job['alignment'])==job['alignment_sha256']
        records=list(SeqIO.parse(job['alignment'],'fasta'));seq={r.id:str(r.seq).replace('X','-') for r in records}
        text=gzip.open(ck,'rt').read();tree=Phylo.read(StringIO(re.search(r'^ newick: (.+)$',text,re.M).group(1)),'newick')
        for tip in tree.get_terminals():tip.name=records[int(tip.name)].id
        selected=sorted([r for r in nodes if r['job_id']==job['job_id']],key=lambda r:int(r['level']));assert len(selected)==3
        sigs=signatures(tree);by={hashlib.sha256(json.dumps(sig,separators=(',',':')).encode()).hexdigest():n for n,sig in sigs.items()}
        assert len(by)==len(sigs);targets=[by[r['partition_signature_sha256']] for r in selected]
        p=params[name];pi=np.array([p['reconstructed_iqtree_empirical_frequencies'][a] for a in AA])
        got,ll=independent(tree,seq,pi,models[job['model'].split('+')[0]],p['gamma_category_rates'],targets)
        with np.load(source/(name+'.npz'),allow_pickle=False) as saved:
            assert np.array_equal(saved['levels'],[0,1,2]) and ''.join(saved['amino_acids'])==AA
            assert got.shape==saved['posterior'].shape and np.isfinite(got).all()
            err=float(np.max(abs(got-saved['posterior'])));le=float(np.max(abs(ll[None,:]-saved['site_log_likelihood'])))
            assert err<plan['posterior_absolute_tolerance'] and le<plan['site_log_likelihood_tolerance'],(name,err,le)
            rows.append(dict(job_id=name,nodes=3,sites=got.shape[1]*3,probabilities=got.size,maximum_probability_difference=err,maximum_site_log_likelihood_difference=le))
        print(name,err,le,flush=True)
    assert len(rows)==78 and sum(r['sites'] for r in rows)==sr['node_sites']
    with (out/'posterior_readback.tsv').open('w') as h:
        w=csv.DictWriter(h,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    verify();receipt=dict(status='passed_all_78_whole_protein_posteriors',fits=78,nodes=234,sites=sr['node_sites'],probabilities=sum(r['probabilities'] for r in rows),maximum_probability_difference=max(r['maximum_probability_difference'] for r in rows),maximum_site_log_likelihood_difference=max(r['maximum_site_log_likelihood_difference'] for r in rows),plan_sha256=sha(pp),source_receipt_sha256=sha(source/'receipt.json'),analytic_two_internal_node_enumeration='passed',artifacts={'posterior_readback.tsv':sha(out/'posterior_readback.tsv')},scope='Separate fixed-root inside/outside algorithm, full columns and direct matrix exponential. Shared empirical matrices, saved parameters, alignment and tree inputs. Tests numerical implementation, not optimization convergence, model adequacy, indels or historical node support.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
