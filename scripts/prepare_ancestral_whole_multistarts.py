#!/usr/bin/env python3
"""Prepare reproducible alternate gamma/branch starts for every whole fit."""
import json
from pathlib import Path
from Bio import Phylo
from prepare_case_ancestral_neighborhoods import sha

def main():
    baseline=Path('results/ancestral/whole-protein-model-fits-20260927-v1')
    old=json.loads(Path('metadata/ancestral_whole_model_fit_plan_20260927.json').read_text())
    source=json.loads((baseline/'receipt.json').read_text())
    inputs=Path('results/ancestral/whole-multistart-inputs-20260927-v1');inputs.mkdir(exist_ok=False)
    audit=Path('results/ancestral/whole-protein-model-readback-20260927-v1/receipt.json')
    ar=json.loads(audit.read_text());assert ar['status']=='complete_78_whole_protein_reports_and_independent_likelihoods'
    assert ar['source_receipt_sha256']==sha(baseline/'receipt.json')
    jobs=[];pins={str(baseline/'receipt.json'):sha(baseline/'receipt.json'),str(audit):sha(audit)}
    for result in source['results']:
        job=result['job'];treepath=baseline/job['job_id']/'fit.treefile'
        assert sha(treepath)==result['artifacts']['fit.treefile'];pins[str(treepath)]=sha(treepath)
        for i,(alpha,scale) in enumerate([(.05,.5),(.5,1.),(2.,2.)]):
            tree=Phylo.read(treepath,'newick');tree.root.branch_length=0.
            for node in tree.find_clades():
                if node is not tree.root:node.branch_length*=scale
            path=inputs/(job['job_id']+'-s'+str(i)+'.nwk')
            Phylo.write(tree,path,'newick',format_branch_length='%1.17g');pins[str(path)]=sha(path)
            restored=Phylo.read(path,'newick')
            assert [(n.name,n.branch_length) for n in restored.find_clades()]==[(n.name,n.branch_length) for n in tree.find_clades()]
            for lower in [.02,.005]:
                j=dict(job);j.update(base_job_id=job['job_id'],job_id=job['job_id']+'-s'+str(i)+'-amin'+str(lower),tree=str(path),tree_sha256=sha(path),start_alpha=alpha,branch_scale=scale,alpha_min=lower,seed=20260927+i);jobs.append(j)
    for f in ['scripts/run_ancestral_whole_multistarts.py','scripts/prepare_ancestral_whole_multistarts.py','scripts/prepare_case_ancestral_neighborhoods.py','data/software_audits/iqtree-3.0.1/rategamma.cpp','data/software_audits/iqtree-3.0.1/tools.cpp',old['iqtree']]:pins[f]=sha(f)
    assert len(jobs)==468
    p=dict(iqtree=old['iqtree'],output='results/ancestral/whole-multistarts-20260927-v1',jobs=jobs,pins=pins,resources=dict(concurrent_fits=2,threads_per_fit=4,aggregate_cpus=8,memory_gib=12,swap_gib=0,output_gib=8,planning_hours=[1,72]),scope='All78 fits x3 initial gamma/branch combinations x2 alpha bounds. Full branch/gamma refitting with epsilon1e-6. No GPU. No convergence or updated posterior claim until independently audited.',gamma_start_semantics_source='https://github.com/iqtree/iqtree3/blob/v3.0.1/model/rategamma.cpp')
    Path('metadata/ancestral_whole_multistart_plan_20260927.json').write_text(json.dumps(p,indent=2)+'\n');print('Prepared',len(jobs),'jobs')

if __name__=='__main__':main()
