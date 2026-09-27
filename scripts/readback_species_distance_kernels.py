#!/usr/bin/env python3
"""Check every tree distance using DendroPy before and after an explicit reroot."""
import argparse,json,hashlib
from pathlib import Path
import dendropy,numpy as np
from scipy.linalg import eigvalsh

sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def distances(tree,taxa):
    pdm=tree.phylogenetic_distance_matrix();lookup={t.label:t for t in tree.taxon_namespace};assert set(lookup)==set(taxa)
    d=np.zeros((len(taxa),len(taxa)))
    for i,a in enumerate(taxa):
        for j in range(i):d[i,j]=d[j,i]=pdm.distance(lookup[a],lookup[taxa[j]])
    return d

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    for p,h in plan['pins'].items():assert sha(p)==h,p
    root=Path(plan['output']);r=json.loads((root/'receipt.json').read_text());assert r['status']=='complete_species_distance_kernels_pending_independent_readback' and r['plan_sha256']==ph
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    taxa=json.loads((root/'taxa.json').read_text());assert len(taxa)==len(set(taxa))==526
    h=np.eye(len(taxa))-np.ones((len(taxa),len(taxa)))/len(taxa);summaries={}
    for label,spec in plan['trees'].items():
        t=dendropy.Tree.get(path=spec['tree'],schema='newick',preserve_underscores=True);d=distances(t,taxa)
        saved=np.load(root/(label+'.npz'));assert set(saved.files)=={'distances','centered_kernel','eigenvalues'}
        assert np.allclose(d,saved['distances'],rtol=1e-11,atol=1e-11)
        kernel=-.5*h@d@h;assert np.allclose(kernel,saved['centered_kernel'],rtol=1e-10,atol=1e-11)
        eigen=eigvalsh(kernel,driver='evr');assert np.allclose(eigen,saved['eigenvalues'],rtol=1e-9,atol=1e-10)
        assert eigen.min()>-1e-10*max(1,eigen.max())
        recovered=np.diag(kernel)[:,None]+np.diag(kernel)[None,:]-2*kernel;assert np.allclose(recovered,d,rtol=1e-10,atol=1e-10)
        edge=t.find_node_with_taxon_label(taxa[-1]).edge;length=float(edge.length)
        t.reroot_at_edge(edge,length1=length/2,length2=length/2,suppress_unifurcations=True)
        rerooted=distances(t,taxa);assert np.allclose(d,rerooted,rtol=1e-11,atol=1e-11)
        assert np.allclose(-.5*h@rerooted@h,kernel,rtol=1e-10,atol=1e-10)
        summaries[label]=dict(unordered_tip_pairs=len(taxa)*(len(taxa)-1)//2,maximum_distance_error=float(abs(d-saved['distances']).max()),maximum_kernel_error=float(abs(kernel-saved['centered_kernel']).max()),maximum_reroot_distance_error=float(abs(d-rerooted).max()),minimum_kernel_eigenvalue=float(eigen.min()),reroot_tip=taxa[-1])
        print('Verified all distances and reroot',label,flush=True)
    for p,hh in plan['pins'].items():assert sha(p)==hh,p
    result=dict(status='passed_full_species_distance_kernel_readback',taxa=len(taxa),trees=summaries,tip_pairs_checked=sum(x['unordered_tip_pairs'] for x in summaries.values()),source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),scope='Every unordered tip-pair distance independently reconstructed with DendroPy and again after explicit terminal-edge rerooting. All centered-kernel entries, recovered distances and eigenvalues checked with alternative linear algebra. Mathematical geometry only, not estimated evolutionary covariance, time calibration or phylogenetic model adequacy.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
