#!/usr/bin/env python3
"""Select best audited feasible starts, including baseline, at each gamma bound."""
import csv,json
from pathlib import Path
from prepare_case_ancestral_neighborhoods import sha

def main():
    old=json.loads(Path('metadata/ancestral_domain_multistart_plan_20260927.json').read_text())
    roots=[Path('results/ancestral/case-domain-model-fits-20260927-v1'),Path(old['output'])]
    audits=[Path('results/ancestral/case-domain-model-readback-20260927-v1'),Path('results/ancestral/domain-multistart-readback-20260927-v1')]
    pins={};candidates={};jobs=[]
    for root,audit in zip(roots,audits):
        ar=json.loads((audit/'receipt.json').read_text());assert ar['source_receipt_sha256']==sha(root/'receipt.json')
        source=json.loads((root/'receipt.json').read_text());results={r['job']['job_id']:r for r in source['results']}
        for f in [root/'receipt.json',audit/'receipt.json',audit/'fit_readback.tsv',audit/'parameters.json']:pins[str(f)]=sha(f)
        for name in ['fit_readback.tsv','parameters.json']:assert sha(audit/name)==ar['artifacts'][name]
        parameters={r['job_id']:r for r in json.loads((audit/'parameters.json').read_text())}
        for row in csv.DictReader((audit/'fit_readback.tsv').open(),delimiter='\t'):
            r=results[row['job_id']];job=r['job'];base=job.get('base_job_id',job['job_id']);alpha=parameters[row['job_id']]['gamma_shape'];ll=float(row['checkpoint_log_likelihood'])
            for lower in [.02,.005]:
                if alpha<lower or ('alpha_min' in job and job['alpha_min']!=lower):continue
                candidates.setdefault((base,lower),[]).append((ll,row['job_id'],root,r,alpha))
    assert len(candidates)==312
    for (base,lower),values in sorted(candidates.items()):
        ll,name,root,result,alpha=max(values,key=lambda x:(x[0],x[1]));tree=root/name/'fit.treefile';assert sha(tree)==result['artifacts']['fit.treefile'];pins[str(tree)]=sha(tree)
        j=dict(result['job']);j.update(job_id=base+'-refine-amin'+str(lower),base_job_id=base,tree=str(tree),tree_sha256=sha(tree),start_alpha=alpha,branch_scale=1.,alpha_min=lower,seed=20260928,selected_start_job_id=name,selected_start_root=str(root),selected_start_log_likelihood=ll,eligible_audited_starts=len(values));jobs.append(j)
    for f in ['scripts/prepare_ancestral_domain_refinements.py','scripts/run_ancestral_domain_refinements.py','scripts/prepare_case_ancestral_neighborhoods.py','metadata/ancestral_domain_multistart_audit_completed_20260927.json','metadata/ancestral_domain_likelihood_replay_completed_20260927.json',old['iqtree']]:pins[f]=sha(f)
    p=dict(iqtree=old['iqtree'],output='results/ancestral/domain-refinements-20260927-v1',jobs=jobs,pins=pins,resources=dict(concurrent_fits=2,threads_per_fit=4,aggregate_cpus=8,memory_gib=12,swap_gib=0,output_gib=4,planning_hours=[.5,24]),scope='312 best feasible starts from audited baseline and same-bound multistarts. Tighter epsilon1e-8; full branch/gamma refitting. No global convergence claim. Numerical audit and posterior update required.')
    Path('metadata/ancestral_domain_refinement_plan_20260927.json').write_text(json.dumps(p,indent=2)+'\n');print('Prepared',len(jobs),'refinements; baseline starts',sum(j['selected_start_root']==str(roots[0]) for j in jobs))

if __name__=='__main__':main()
