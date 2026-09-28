#!/usr/bin/env python3
"""Whole-protein amino-acid marginals at topology-matched candidate ancestors."""
import csv
import subprocess
import time
import psutil
import gzip
import hashlib
import json
import re
from io import StringIO
from pathlib import Path
import numpy as np
from Bio import Phylo,SeqIO
from scipy.special import xlogy
from threadpoolctl import threadpool_limits
from infer_refined_domain_ancestors import posterior
from replay_ancestral_domain_likelihoods import AA,matrices,generator
from map_ancestral_fitted_nodes import signatures
from prepare_case_ancestral_neighborhoods import sha


def main():
    pp=Path('metadata/refined_whole_ancestor_plan_v2_20260927.json');plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():
        assert sha(p)==h,p
    identity=plan['auditor']
    while True:
        state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',identity['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus','-p','MainPID'],text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:
            assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0',state
            break
        assert state['ActiveState'] in ['active','activating','deactivating'],state
        try:
            process=psutil.Process(identity['pid'])
            assert process.create_time()==identity['created'] and process.cmdline()==identity['cmdline'] and int(state['MainPID'])==process.pid
        except psutil.NoSuchProcess:
            time.sleep(1);continue
        time.sleep(30)
    for p,h in plan['pins'].items():assert sha(p)==h,p
    fits=Path(plan['fits']);audit_path=Path(plan['audit_receipt'])
    audit=json.loads(audit_path.read_text())
    assert audit['status']=='complete_156_refinement_reports_and_independent_likelihoods'
    assert audit['fits']==156 and audit['audit_plan_sha256']==sha(plan['audit_plan'])
    assert audit['source_receipt_sha256']==sha(fits/'receipt.json')
    for name,h in audit['artifacts'].items():
        assert sha(audit_path.parent/name)==h
    assert sha(plan['parameters'])==audit['artifacts']['parameters.json']
    params={r['job_id']:r for r in json.loads(Path(plan['parameters']).read_text())}
    source=json.loads((fits/'receipt.json').read_text())
    assert len(source['results'])==156 and len(params)==156
    assert {f['job']['job_id'] for f in source['results']}==set(params)
    with Path(plan['source_nodes']).open() as handle:
        source_nodes=[r for r in csv.DictReader(handle,delimiter='\t') if r['guide']=='profile' and r['dataset']=='whole']
    baseline_jobs={j['job_id']:j for j in json.loads(Path(plan['baseline_plan']).read_text())['jobs']}
    models=matrices(plan['exchangeabilities']);out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    summaries=[];mappings=[];sites=0
    for fit in source['results']:
        job=fit['job'];name=job['job_id'];folder=fits/name
        assert sha(folder/'fit.ckp.gz')==fit['artifacts']['fit.ckp.gz']
        assert sha(job['alignment'])==job['alignment_sha256']
        assert sha(job['tree'])==job['tree_sha256']
        records=list(SeqIO.parse(job['alignment'],'fasta'))
        original={r.id:str(r.seq) for r in records};assert len(original)==len(records)
        assert set(''.join(original.values()))<=set(AA+'-X')
        sequences={g:s.replace('X','-') for g,s in original.items()}
        text=gzip.open(folder/'fit.ckp.gz','rt').read()
        tree=Phylo.read(StringIO(re.search(r'^ newick: (.+)$',text,re.M).group(1)),'newick')
        for tip in tree.get_terminals():
            tip.name=records[int(tip.name)].id
        expected=float(re.search(r'^ 0: ([-+0-9.eE]+)',text,re.M).group(1))
        param=params[name];pi=np.array([param['reconstructed_iqtree_empirical_frequencies'][a] for a in AA])
        q=generator(models[job['model'].split('+')[0]],pi);rates=param['gamma_category_rates']
        baseline=baseline_jobs[job['base_job_id']];assert sha(baseline['tree'])==baseline['tree_sha256']
        guide=Phylo.read(baseline['tree'],'newick');guide_sigs=signatures(guide);fitted_sigs=signatures(tree)
        assert {n.name for n in guide.get_terminals()}==set(sequences)
        arrays=[];likelihoods=[]
        for row in sorted([r for r in source_nodes if r['family']==job['family']],key=lambda r:int(r['level'])):
            level=int(row['level']);names=json.loads(row['output_nodes_json']);assert len(names)==1
            matches=[n for n in guide_sigs if n.name==names[0]];assert len(matches)==1
            sig=guide_sigs[matches[0]];targets=[n for n in fitted_sigs if fitted_sigs[n]==sig]
            signature=hashlib.sha256(json.dumps(sig,separators=(',',':')).encode()).hexdigest()
            if len(sig)==2:
                assert level==3 and not targets
                mappings.append(dict(job_id=name,level=level,status='degree_two_root_position_not_identified',partition_signature_sha256=signature))
                continue
            assert len(targets)==1 and level in [0,1,2]
            target=targets[0];prob,logs=posterior(tree,target,sequences,pi,q,rates)
            assert prob.shape==(job['columns'],20) and np.isfinite(prob).all()
            error=float(logs.sum()-expected);assert abs(error)<.001,(name,level,error)
            arrays.append(prob);likelihoods.append(logs);sites+=job['columns']
            mappings.append(dict(job_id=name,level=level,status='matched_unrooted_vertex',partition_signature_sha256=signature))
            summaries.append(dict(job_id=name,level=level,columns=job['columns'],likelihood_difference=error,
                mean_maximum_probability=float(prob.max(axis=1).mean()),sites_below_090=int((prob.max(axis=1)<.9).sum()),
                mean_entropy_nats=float((-xlogy(prob,prob).sum(axis=1)).mean()),
                unknown_X_input_cells=sum(s.count('X') for s in original.values())))
        assert len(arrays)==3 and max(np.max(abs(x-likelihoods[0])) for x in likelihoods)<1e-8
        np.savez_compressed(out/(name+'.npz'),posterior=np.array(arrays),site_log_likelihood=np.array(likelihoods),levels=np.array([0,1,2]),amino_acids=np.array(list(AA)))
        print(name,'three_whole_protein_ancestors_produced',flush=True)
    assert len(summaries)==468 and len(mappings)==624
    for name,rows in [('node_summary.tsv',summaries),('candidate_mapping.tsv',mappings)]:
        with (out/name).open('w') as handle:
            writer=csv.DictWriter(handle,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    for p,h in plan['pins'].items():
        assert sha(p)==h,p
    receipt=dict(status='all_156_refined_whole_protein_marginals_produced_pending_independent_audit',fits=156,nodes=468,
        node_sites=sites,probabilities=sites*20,plan_sha256=sha(pp),audit_receipt_sha256=sha(audit_path),
        maximum_likelihood_difference=max(abs(r['likelihood_difference']) for r in summaries),
        artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},scope=plan['scope'])
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    with threadpool_limits(limits=1):
        main()
