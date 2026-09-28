#!/usr/bin/env python3
"""Refine best audited whole-protein fits under both gamma lower bounds."""
import argparse,concurrent.futures,fcntl,json,math,re,subprocess,time
from pathlib import Path
from Bio import Phylo,SeqIO
from prepare_case_ancestral_neighborhoods import sha


def splits(tree):
    alltips={n.name for n in tree.get_terminals()};desc={};result=set()
    for node in tree.find_clades(order='postorder'):
        side=set().union(*(desc[c] for c in node.clades)) if node.clades else {node.name}
        desc[node]=side;other=alltips-side
        if min(len(side),len(other))>1:
            result.add(min((tuple(sorted(side)),tuple(sorted(other))),key=lambda x:(len(x),x)))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());digest=sha(a.plan)
    def verify():
        assert sha(a.plan)==digest
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True)
    lock=(out/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    cp=out/'plan.json'
    if cp.exists():assert json.loads(cp.read_text())==plan
    else:cp.write_text(json.dumps(plan,indent=2)+'\n')
    jobs=plan['jobs'];assert len(jobs)==156 and len({r['job_id'] for r in jobs})==156
    def run(job):
        folder=out/job['job_id'];rp=folder/'receipt.json'
        if rp.exists():
            r=json.loads(rp.read_text());assert r['status']=='complete_refined_whole_fit_pending_full_audit' and r['plan_sha256']==digest
            for name,h in r['artifacts'].items():assert sha(folder/name)==h
            return r
        folder.mkdir(exist_ok=False)
        alignment=Path(job['alignment']);topology=Path(job['tree'])
        assert sha(alignment)==job['alignment_sha256'] and sha(topology)==job['tree_sha256']
        records=list(SeqIO.parse(alignment,'fasta'));tips={r.id for r in records}
        assert len(records)==len(tips)==job['proteins'] and {len(r.seq) for r in records}=={job['columns']}
        source=Phylo.read(topology,'newick');assert {n.name for n in source.get_terminals()}==tips
        command=[plan['iqtree'],'-s',str(alignment.resolve()),'-st','AA','-t',str(topology.resolve()),'--tree-fix','-m',job['model'],'-T','4','--mem','4G','--seed',str(job['seed']),'-a',str(job['start_alpha']),'-optfromgiven','--alpha-min',str(job['alpha_min']),'--epsilon','0.00000001','-keep-ident','--prefix',str((folder/'fit').resolve())]
        started=time.monotonic()
        with (folder/'stdout.log').open('a') as handle:subprocess.run(command,stdout=handle,stderr=subprocess.STDOUT,check=True)
        fitted=Phylo.read(folder/'fit.treefile','newick');ftips=[n.name for n in fitted.get_terminals()]
        assert set(ftips)==tips and len(ftips)==len(tips) and splits(fitted)==splits(source)
        for node in fitted.find_clades():
            if node is not fitted.root:assert node.branch_length is not None and math.isfinite(node.branch_length) and node.branch_length>=0
        report=(folder/'fit.iqtree').read_text()
        likelihood=re.search(r'^Log-likelihood of the tree:\s+([-+0-9.eE]+)',report,re.M)
        assert likelihood and math.isfinite(float(likelihood.group(1)))
        assert 'Model of substitution: '+job['model'] in report
        r=dict(status='complete_refined_whole_fit_pending_full_audit',job=job,plan_sha256=digest,command=command,elapsed_seconds=time.monotonic()-started,reported_log_likelihood=float(likelihood.group(1)),tips=len(tips),unrooted_internal_splits=len(splits(fitted)),artifacts={p.name:sha(p) for p in folder.iterdir()})
        rp.write_text(json.dumps(r,indent=2)+'\n');print(job['job_id'],'complete',r['reported_log_likelihood'],flush=True);return r
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,jobs))
    verify()
    result=dict(status='complete_156_best_start_whole_refinements_pending_independent_audit',plan_sha256=digest,results=results,scope='Best audited feasible solution per78 baseline models and two gamma bounds, restarted with epsilon1e-8 and free branch/gamma parameters. Completion alone does not establish convergence or qualify updated ancestral probabilities.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
