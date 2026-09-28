#!/usr/bin/env python3
"""Read back all fixed-topology model reports, checkpoint parameters and warnings."""
import csv,gzip,json,math,re,subprocess,time
import psutil
from threadpoolctl import threadpool_limits
from replay_ancestral_domain_likelihoods import AA,matrices,generator,likelihood as replay_likelihood
from collections import Counter,defaultdict
from io import StringIO
from pathlib import Path
import numpy as np
from scipy.special import gammainc
from scipy.stats import gamma
from Bio import Phylo,SeqIO
from prepare_case_ancestral_neighborhoods import sha


def edge_map(tree):
    allnames={t.name for t in tree.get_terminals()};desc={};edges={}
    for node in tree.find_clades(order='postorder'):
        side=set().union(*(desc[c] for c in node.clades)) if node.clades else {node.name}
        desc[node]=side
        if node is tree.root:continue
        other=allnames-side
        key=min((tuple(sorted(side)),tuple(sorted(other))),key=lambda x:(len(x),x))
        assert node.branch_length is not None and math.isfinite(node.branch_length) and node.branch_length>=0
        edges[key]=edges.get(key,0.)+node.branch_length
    return edges


def main():
    audit_pp=Path('metadata/ancestral_whole_refinement_audit_plan_20260927.json');audit_plan=json.loads(audit_pp.read_text())
    def verify_audit():
        for p,h in audit_plan['pins'].items():assert sha(p)==h,p
    verify_audit();identity=audit_plan['producer']
    while True:
        state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',identity['unit'],'--property=ActiveState,Result,ExecMainStatus,MainPID'],text=True).splitlines())
        if state['ActiveState']=='inactive':
            assert state['Result']=='success' and state['ExecMainStatus']=='0';break
        assert state['ActiveState'] in ['active','activating','deactivating'],state
        proc=psutil.Process(identity['pid']);assert proc.create_time()==identity['created'] and proc.cmdline()==identity['cmdline'] and int(state['MainPID'])==proc.pid
        print('waiting_for_verified_refinement_producer',flush=True);time.sleep(30)
    verify_audit();models=matrices(audit_plan['exchangeabilities'])
    pp=Path('metadata/ancestral_whole_refinement_plan_20260927.json');plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():assert sha(p)==h
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',identity['unit'],'--property=ActiveState,Result,ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    root=Path(plan['output']);rp=root/'receipt.json';source=json.loads(rp.read_text())
    assert source['status']=='complete_156_best_start_whole_refinements_pending_independent_audit' and source['plan_sha256']==sha(pp)
    runs={r['job']['job_id']:r for r in source['results']};assert len(runs)==len(plan['jobs'])==156
    rows=[];warnings=[];parameters=[]
    for job in plan['jobs']:
        name=job['job_id'];folder=root/name;receipt=json.loads((folder/'receipt.json').read_text());assert receipt==runs[name] and receipt['job']==job
        for p,h in receipt['artifacts'].items():assert sha(folder/p)==h
        assert sha(job['alignment'])==job['alignment_sha256'] and sha(job['tree'])==job['tree_sha256']
        records=list(SeqIO.parse(job['alignment'],'fasta'));tips={r.id for r in records}
        assert len(records)==len(tips)==job['proteins'] and {len(r.seq) for r in records}=={job['columns']}
        report=(folder/'fit.iqtree').read_text();log=(folder/'fit.log').read_text();checkpoint=gzip.open(folder/'fit.ckp.gz','rt').read()
        assert re.search(r'^finished: true$',checkpoint,re.M)
        assert 'Model of substitution: '+job['model'] in report
        def numeric(pattern,text=report):
            match=re.search(pattern,text,re.M);assert match,(name,pattern);value=float(match.group(1));assert math.isfinite(value);return value
        likelihood=numeric(r'^Log-likelihood of the tree:\s+([-+0-9.eE]+)')
        checkpoint_lh=numeric(r'^ 0:\s+([-+0-9.eE]+)',checkpoint)
        assert abs(likelihood-checkpoint_lh)<5.1e-5 and likelihood==receipt['reported_log_likelihood']
        alpha=numeric(r'^ gamma_shape:\s+([-+0-9.eE]+)',checkpoint);assert alpha>0
        assert abs(alpha-numeric(r'^Gamma shape alpha:\s+([-+0-9.eE]+)'))<5.1e-5
        frequencies=dict((aa,float(v)) for aa,v in re.findall(r'pi\(([A-Z])\) = ([-+0-9.eE]+)',report));assert len(frequencies)==20
        counts=Counter(''.join(str(r.seq) for r in records));canonical='ARNDCQEGHILKMFPSTWYV';total=sum(counts[a] for a in canonical)
        assert set(counts)<=set(canonical+'-X')
        # IQ-TREE 3.0.1 uses eight ambiguity-allocation iterations from
        # uniform frequencies, including gaps; keep_zero_freq defaults true.
        cells=len(records)*job['columns'];estimated={a:1/20 for a in canonical}
        for iteration in range(8):estimated={a:(counts[a]+(counts['-']+counts['X'])*estimated[a])/cells for a in canonical}
        frequency_error=max(abs(frequencies[a]-estimated[a]) for a in canonical);assert frequency_error<5.2e-5,(name,frequency_error)
        rates=re.findall(r'^\s+([1-4])\s+([0-9.]+)\s+(0\.2500)\s*$',report,re.M);assert len(rates)==4
        bounds=gamma.ppf(np.linspace(0,1,5),a=alpha,scale=1/alpha)
        expected=4*np.diff(gammainc(alpha+1,alpha*bounds))
        rate_error=max(abs(float(row[1])-value) for row,value in zip(rates,expected));assert rate_error<5.2e-5
        fitted=Phylo.read(folder/'fit.treefile','newick');ftips=[t.name for t in fitted.get_terminals()];assert len(ftips)==len(tips) and set(ftips)==tips
        edges=edge_map(fitted);original=edge_map(Phylo.read(job['tree'],'newick'));assert set(edges)==set(original)
        checkpoint_tree=re.search(r'^ newick: (.+)$',checkpoint,re.M);assert checkpoint_tree
        ct=Phylo.read(StringIO(checkpoint_tree.group(1)),'newick')
        for tip in ct.get_terminals():tip.name=records[int(tip.name)].id
        ce=edge_map(ct);assert set(ce)==set(edges)
        branch_error=max(abs(ce[k]-edges[k]) for k in edges);assert branch_error<1.1e-9
        k=numeric(r'^Number of free parameters .*?:\s+(\d+)');assert k==len(edges)+20,(name,k,len(edges))
        aic=numeric(r'^Akaike information criterion \(AIC\) score:\s+([-+0-9.eE]+)')
        bic=numeric(r'^Bayesian information criterion \(BIC\) score:\s+([-+0-9.eE]+)')
        assert abs(aic-(2*k-2*checkpoint_lh))<1.1e-4 and abs(bic-(math.log(job['columns'])*k-2*checkpoint_lh))<1.1e-4
        assert abs(sum(edges.values())-numeric(r'^Total tree length .*?:\s+([-+0-9.eE]+)'))<5.2e-5
        for origin,text in [('report',report),('log',log)]:
            for line in text.splitlines():
                if re.search(r'WARNING|ERROR|failed|not converg',line,re.I):warnings.append(dict(job_id=name,origin=origin,message=line))
        parameters.append(dict(job_id=name,gamma_shape=alpha,reconstructed_iqtree_empirical_frequencies=estimated,raw_nongap_count_frequencies={a:counts[a]/total for a in canonical},gamma_category_rates=expected.tolist()))
        pi=np.array([estimated[a] for a in AA]);q=generator(models[job['model'].split('+')[0]],pi)
        computed=replay_likelihood(ct,{r.id:str(r.seq).replace('X','-') for r in records},pi,q,expected)
        replay_error=computed-checkpoint_lh;assert abs(replay_error)<.001,(name,replay_error)
        assert alpha>=job['alpha_min']-1e-10
        assert counts['X']==job['unknown_X_residues']
        selected_folder=Path(job['selected_start_root'])/job['selected_start_job_id']
        selected_receipt=json.loads((selected_folder/'receipt.json').read_text())
        assert sha(selected_folder/'fit.ckp.gz')==selected_receipt['artifacts']['fit.ckp.gz']
        selected_checkpoint=gzip.open(selected_folder/'fit.ckp.gz','rt').read()
        selected_ll=numeric(r'^ 0:\s+([-+0-9.eE]+)',selected_checkpoint)
        assert selected_ll==job['selected_start_log_likelihood']
        rows.append(dict(job_id=name,unknown_X_residues=counts['X'],selected_start_job_id=job['selected_start_job_id'],selected_start_log_likelihood=selected_ll,refinement_likelihood_change=checkpoint_lh-selected_ll,gamma_shape_change=alpha-job['start_alpha'],maximum_branch_change=max(abs(edges[e]-original[e]) for e in edges),base_job_id=job['base_job_id'],start_alpha=job['start_alpha'],branch_scale=job['branch_scale'],alpha_min=job['alpha_min'],recomputed_log_likelihood=computed,likelihood_replay_error=replay_error,family=job['family'],boundary=job['boundary'],method=job['method'],model=job['model'],proteins=len(tips),columns=job['columns'],checkpoint_log_likelihood=checkpoint_lh,free_parameters=int(k),reported_aic=aic,reported_bic=bic,gamma_shape=alpha,frequency_rounding_error=frequency_error,rate_rounding_error=rate_error,branch_serialization_error=branch_error,edges=len(edges),edges_below_1e_5=sum(v<1e-5 for v in edges.values()),edges_above_5=sum(v>5 for v in edges.values()),warning_lines=sum(w['job_id']==name for w in warnings)))
        print(name,'full_report_and_likelihood_audit_passed',flush=True)
    out=Path(audit_plan['output']);out.mkdir(exist_ok=False)
    verify_audit()
    for filename,data in [('fit_readback.tsv',rows),('warnings.tsv',warnings)]:
        with (out/filename).open('w') as f:w=csv.DictWriter(f,list(data[0]) if data else ['job_id','origin','message'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    (out/'parameters.json').write_text(json.dumps(parameters,indent=2)+'\n')
    result=dict(status='complete_156_refinement_reports_and_independent_likelihoods',maximum_likelihood_replay_error=max(abs(r['likelihood_replay_error']) for r in rows),fits_with_likelihood_change_above_001=sum(r['refinement_likelihood_change']>.001 for r in rows),fits_with_likelihood_change_below_minus001=sum(r['refinement_likelihood_change']<-.001 for r in rows),maximum_likelihood_change=max(r['refinement_likelihood_change'] for r in rows),minimum_likelihood_change=min(r['refinement_likelihood_change'] for r in rows),audit_plan_sha256=sha(audit_pp),fits=len(rows),warning_lines=len(warnings),fits_with_warning_lines=sum(r['warning_lines']>0 for r in rows),fits_with_edges_below_1e_5=sum(r['edges_below_1e_5']>0 for r in rows),source_receipt_path=str(rp),source_receipt_sha256=sha(rp),plan_sha256=sha(pp),frequency_implementation_source=dict(url='https://github.com/iqtree/iqtree3/blob/v3.0.1/alignment/alignment.cpp',local_path='data/software_audits/iqtree-3.0.1/alignment.cpp',sha256=sha('data/software_audits/iqtree-3.0.1/alignment.cpp')),script_sha256=sha(__file__),terminal_state=state,artifacts={p.name:sha(p) for p in out.iterdir()},scope='All 156 refinement model labels, tip sets, unrooted edge sets, checkpoint/report likelihood agreement, checkpoint/tree branch lengths, version-specific eight-iteration empirical frequencies, gamma-category means and AIC/BIC arithmetic checked. Warnings and short edges retained. Independent scaled-pruning likelihood replay checked for all156; signed likelihood, gamma and branch changes from selected starts retained. Agreement is not a proof of global convergence or model adequacy; ancestral probabilities remain from the original fits until separately propagated.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
