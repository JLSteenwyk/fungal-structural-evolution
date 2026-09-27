#!/usr/bin/env python3
"""Compute conditional node marginals while retaining alignment/gap uncertainty."""
import csv,gzip,hashlib,json,math,re
from io import StringIO
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,expm
from scipy.special import logsumexp,xlogy
from Bio import Phylo,SeqIO
from threadpoolctl import threadpool_limits
from prepare_case_ancestral_neighborhoods import sha,read
from replay_ancestral_domain_likelihoods import AA,matrices,generator
from map_ancestral_fitted_nodes import signatures


def posterior(tree,target,sequences,pi,q,rates):
    names=list(sequences);array=np.array([list(sequences[g]) for g in names]).T
    patterns,inverse=np.unique(array,axis=0,return_inverse=True)
    evidence={g:np.array([[float(x=='-' or x==a) for x in patterns[:,i]] for a in AA]) for i,g in enumerate(names)}
    graph={n:[] for n in tree.find_clades()}
    for n in tree.find_clades():
        for c in n.clades:graph[n].append((c,c.branch_length));graph[c].append((n,c.branch_length))
    parent={target:None};order=[target]
    for node in order:
        for neighbor,length in graph[node]:
            if neighbor is parent[node]:continue
            assert neighbor not in parent;parent[neighbor]=node;order.append(neighbor)
    assert len(order)==len(graph)
    rootpi=np.sqrt(pi);vals,vecs=eigh(rootpi[:,None]*q/rootpi[None,:]);joints=[]
    for rate in rates:
        values={};scales={}
        for node in reversed(order):
            part=evidence[node.name].copy() if node.name in evidence else np.ones((20,len(patterns)))
            scale=np.zeros(len(patterns))
            for child,length in graph[node]:
                if child is parent[node]:continue
                t=length*rate
                p=expm(q*t) if t<1e-8 else ((vecs*np.exp(vals*t))@vecs.T)*rootpi[None,:]/rootpi[:,None]
                assert p.min()>-1e-12 and np.max(abs(p.sum(axis=1)-1))<1e-10;p=np.maximum(p,0.)
                part*=p@values.pop(child);scale+=scales.pop(child)
                norm=part.max(axis=0);assert (norm>0).all();part/=norm;scale+=np.log(norm)
            values[node]=part;scales[node]=scale
        with np.errstate(divide='ignore'):joint=np.log(pi[:,None]*values[target])+scales[target][None,:]-math.log(len(rates))
        joints.append(joint)
    state_logs=logsumexp(np.array(joints),axis=0);site=logsumexp(state_logs,axis=0)
    marginal=np.exp(state_logs-site[None,:]);assert np.max(abs(marginal.sum(axis=0)-1))<1e-10
    return marginal[:,inverse].T,site[inverse]


def main():
    pp=Path('metadata/conditional_domain_ancestor_plan_20260927.json');plan=json.loads(pp.read_text())
    def verify():
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();models=matrices(plan['exchangeabilities'])
    pi=np.arange(1,21,dtype=float);pi/=pi.sum();q=generator(models['LG'],pi)
    toy=Phylo.read(StringIO('(a:0.07,b:0.13,c:0.11)n0;'),'newick');seq={'a':'A-','b':'NW','c':'R-'};rates=[.2,.6,1.1,2.1]
    observed,site=posterior(toy,toy.root,seq,pi,q,rates)
    exact=[]
    for col in range(2):
        joint=np.zeros(20)
        for rate in rates:
            v=pi.copy()
            for node in toy.get_terminals():
                state=seq[node.name][col]
                if state!='-':v*=expm(q*node.branch_length*rate)[:,AA.index(state)]
            joint+=v/4
        exact.append(joint/joint.sum());assert abs(site[col]-math.log(joint.sum()))<1e-11
    assert np.max(abs(observed-np.array(exact)))<1e-11
    unknown,_=posterior(toy,toy.root,{g:'-' for g in seq},pi,q,rates);assert np.max(abs(unknown[0]-pi))<1e-11
    fits=Path(plan['fits']);source=json.loads((fits/'receipt.json').read_text())
    parameters={p['job_id']:p for p in json.loads(Path(plan['parameters']).read_text())}
    nodes=read(plan['nodes']);nodes=[n for n in nodes if n['guide']=='profile' and n['status']=='matched_unrooted_vertex'];assert len(nodes)==468
    source_nodes={(r['family'],int(r['level'])):set(json.loads(r['retained_set_json'])) for r in read(plan['source_nodes']) if r['guide']=='profile' and r['dataset']=='domain'}
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);summary=[];site_rows=0
    with gzip.open(out/'site_marginals.tsv.gz','wt') as handle:
        writer=None
        for fit in source['results']:
            job=fit['job'];name=job['job_id'];folder=fits/name;assert sha(folder/'fit.ckp.gz')==fit['artifacts']['fit.ckp.gz'];assert sha(job['alignment'])==job['alignment_sha256']
            records=list(SeqIO.parse(job['alignment'],'fasta'));sequences={r.id:str(r.seq) for r in records};assert len(sequences)==len(records)
            text=gzip.open(folder/'fit.ckp.gz','rt').read();tree=Phylo.read(StringIO(re.search(r'^ newick: (.+)$',text,re.M).group(1)),'newick')
            for tip in tree.get_terminals():tip.name=records[int(tip.name)].id
            expected_lh=float(re.search(r'^ 0: ([-+0-9.eE]+)',text,re.M).group(1))
            param=parameters[name];pi=np.array([param['reconstructed_iqtree_empirical_frequencies'][a] for a in AA]);q=generator(models[job['model'].split('+')[0]],pi);rates=param['gamma_category_rates']
            sigs=signatures(tree);targets=sorted([n for n in nodes if n['job_id']==name],key=lambda n:int(n['level']));assert len(targets)==3
            arrays=[];loglikes=[]
            for record in targets:
                target=[n for n in sigs if n.name==record['fitted_node']];assert len(target)==1;target=target[0]
                assert hashlib.sha256(json.dumps(sigs[target],separators=(',',':')).encode()).hexdigest()==record['partition_signature_sha256']
                probabilities,logs=posterior(tree,target,sequences,pi,q,rates)
                assert probabilities.shape==(job['columns'],20) and np.isfinite(probabilities).all()
                error=float(logs.sum())-expected_lh;assert abs(error)<.001,(name,record['level'],error)
                arrays.append(probabilities);loglikes.append(logs);level=int(record['level']);desc=source_nodes[job['family'],level];assert desc<=set(sequences)
                best=probabilities.max(axis=1);states=probabilities.argmax(axis=1);entropy=-xlogy(probabilities,probabilities).sum(axis=1)
                for col in range(job['columns']):
                    inside=sum(sequences[g][col]!='-' for g in desc);outside=sum(s[col]!='-' for g,s in sequences.items() if g not in desc)
                    row=dict(job_id=name,family=job['family'],boundary=job['boundary'],method=job['method'],model=job['model'],level=level,fitted_node=target.name,column=col+1,map_amino_acid=AA[states[col]],maximum_probability=float(best[col]),entropy_nats=float(entropy[col]),descendant_observed=inside,descendant_total=len(desc),outside_observed=outside,outside_total=len(sequences)-len(desc))
                    if writer is None:writer=csv.DictWriter(handle,list(row),delimiter='\t',lineterminator='\n');writer.writeheader()
                    writer.writerow(row);site_rows+=1
                summary.append(dict(job_id=name,level=level,fitted_node=target.name,columns=job['columns'],likelihood_difference=error,mean_maximum_probability=float(best.mean()),sites_below_090=int((best<.9).sum()),mean_entropy_nats=float(entropy.mean())))
            assert max(np.max(abs(x-loglikes[0])) for x in loglikes)<1e-8
            path=out/(name+'.npz');np.savez_compressed(path,posterior=np.array(arrays),site_log_likelihood=np.array(loglikes),levels=np.array([0,1,2]),amino_acids=np.array(list(AA)))
            with np.load(path,allow_pickle=False) as saved:assert np.array_equal(saved['posterior'],np.array(arrays)) and np.array_equal(saved['site_log_likelihood'],np.array(loglikes))
            print(name,'three_internal_nodes_complete',flush=True)
    assert len(summary)==468
    with (out/'node_summary.tsv').open('w') as handle:w=csv.DictWriter(handle,list(summary[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summary)
    verify();r=dict(status='complete_conditional_domain_node_marginals_pending_independent_posterior_audit',fits=156,nodes=468,site_rows=site_rows,plan_sha256=sha(pp),analytic_three_tip_and_unknown_checks='passed',maximum_likelihood_difference=max(abs(x['likelihood_difference']) for x in summary),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All 20 amino-acid marginals at three mapped internal vertices for every fitted model/alignment. Conditional on saved parameters/topology; gaps unknown, not ancestral deletions. Complete node/site uncertainty retained. No final ancestral FASTA, ensemble, GPU prediction, model adequacy or optimization convergence claim.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='artifacts'}),flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
