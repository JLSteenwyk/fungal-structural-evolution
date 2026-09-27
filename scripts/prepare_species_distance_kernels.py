#!/usr/bin/env python3
"""Derive full patristic distances and centered tree kernels from audited species trees."""
import argparse,json,math
from pathlib import Path
import numpy as np
from Bio import Phylo
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);taxa=None;summaries={}
    for label,spec in plan['trees'].items():
        receipt=json.loads(Path(spec['receipt']).read_text());audit=json.loads(Path(spec['audit']).read_text());assert audit['status']==spec['audit_status']
        if spec['audit']!=spec['receipt']:assert audit['source_receipt_sha256']==sha(spec['receipt'])
        assert sha(spec['tree'])==receipt['artifacts'][Path(spec['tree']).name]
        tree=Phylo.read(spec['tree'],'newick');tips=tree.get_terminals();names=sorted(t.name for t in tips);assert len(names)==len(set(names))==526
        if taxa is None:taxa=names
        assert taxa==names;index={name:i for i,name in enumerate(taxa)}
        nodes=list(tree.find_clades(order='preorder'));edges=[n for n in nodes if n is not tree.root];lengths=np.array([float(n.branch_length) for n in edges]);assert np.isfinite(lengths).all() and (lengths>=0).all()
        incidence=np.zeros((len(taxa),len(edges)),dtype=float)
        for j,node in enumerate(edges):
            for tip in node.get_terminals():incidence[index[tip.name],j]=1
        features=incidence*np.sqrt(lengths)[None,:];gram=features@features.T;depth=np.diag(gram);dist=depth[:,None]+depth[None,:]-2*gram
        np.fill_diagonal(dist,0);assert dist.min()>-1e-10
        centered=features-features.mean(axis=0);kernel=centered@centered.T
        independent=-.5*(dist-dist.mean(axis=0)[None,:]-dist.mean(axis=1)[:,None]+dist.mean())
        assert np.allclose(kernel,independent,rtol=1e-12,atol=1e-12) and np.max(abs(kernel.sum(axis=0)))<1e-8
        eig=np.linalg.eigvalsh(kernel);assert eig.min()>-1e-10*max(1,eig.max())
        np.savez_compressed(out/(label+'.npz'),distances=dist,centered_kernel=kernel,eigenvalues=eig)
        summaries[label]=dict(taxa=len(taxa),edges=len(edges),total_branch_length=float(lengths.sum()),maximum_patristic_distance=float(dist.max()),minimum_kernel_eigenvalue=float(eig.min()),maximum_kernel_eigenvalue=float(eig.max()),tree_sha256=sha(spec['tree']),audit_sha256=sha(spec['audit']))
        print('Prepared all-tip kernel',label,flush=True)
    (out/'taxa.json').write_text(json.dumps(taxa,indent=2)+'\n');verify()
    result=dict(status='complete_species_distance_kernels_pending_independent_readback',plan_sha256=ph,trees=summaries,artifacts={p.name:sha(p) for p in out.iterdir()},scope='All 526 tips on five audited completed tree alternatives. Patristic distances from edge-path features; centered kernel equals minus half the doubly centered distance matrix (not squared distances). Root-invariant for zero-sum contrasts. Substitution units, not dated branch durations. This mathematical kernel is not an empirically fitted covariance, model adequacy result or final species phylogeny; fourth PMSF and hybrid sensitivities remain pending.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
