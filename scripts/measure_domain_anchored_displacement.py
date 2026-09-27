#!/usr/bin/env python3
"""Fit domain-contained whole-protein correspondences and measure outside-domain residuals."""
import csv,gzip,json
from collections import Counter,defaultdict
from functools import lru_cache
from pathlib import Path
import numpy as np
from assess_domain_alignment_geometry import geometry
from duplication_alignment_numeric_readback import load_pdb
from screen_duplication_domain_alignment_coverage import sha


def transform(x,y,method):
    xc=x.mean(0);yc=y.mean(0);h=(x-xc).T@(y-yc)
    if method=='svd':
        u,s,vt=np.linalg.svd(h);d=np.eye(3);d[-1,-1]=1 if np.linalg.det(u@vt)>=0 else -1;r=u@d@vt
    else:
        z=np.array([h[1,2]-h[2,1],h[2,0]-h[0,2],h[0,1]-h[1,0]])
        k=np.empty((4,4));k[0,0]=np.trace(h);k[0,1:]=z;k[1:,0]=z;k[1:,1:]=h+h.T-np.trace(h)*np.eye(3)
        _,vectors=np.linalg.eigh(k);w=vectors[0,-1];v=vectors[1:,-1]
        cross=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
        r=((w*w-v@v)*np.eye(3)+2*np.outer(v,v)+2*w*cross).T
    assert np.allclose(r.T@r,np.eye(3),atol=1e-12) and abs(np.linalg.det(r)-1)<1e-12
    return r,yc-xc@r


def residual(x,y,r,t):return float(np.sqrt(np.mean(np.sum((x@r+t-y)**2,axis=1))))


def fixture():
    rng=np.random.default_rng(724);x=rng.normal(size=(60,3));q,_=np.linalg.qr(rng.normal(size=(3,3)))
    if np.linalg.det(q)<0:q[:,-1]*=-1
    for shift in [0.,2.]:
        y=x@q+np.array([3.,-2.,7.]);y[30:]+=np.array([shift,0.,0.])
        for method in ['svd','quaternion']:
            r,t=transform(x[:30],y[:30],method)
            assert residual(x[:30],y[:30],r,t)<1e-12
            assert abs(residual(x[30:],y[30:],r,t)-shift)<1e-12
    return dict(rigid_and_outside_translation_cases=4,status='passed')


def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    fixtures=fixture();sources={}
    def checked(root,name):
        rp=root/'receipt.json';p=root/name;r=json.loads(rp.read_text());assert sha(p)==r['artifacts'][name]
        sources[str(rp)]=sha(rp);sources[str(p)]=sha(p);return p
    base=Path('results/structural_comparisons');case=base/'whole-domain-case-dossiers-20260927-v1'
    domain_links=rows(checked(case,'domain_reference_links.tsv'));whole_links=rows(checked(case,'whole_reference_links.tsv'))
    wm=base/'whole-protein-common-residues-20260927-v1';dm=base/'duplication-domain-common-residues-20260927-v1'
    wtriads={r['triad_id']:r for r in rows(checked(wm,'model_triads.tsv'))}
    lookup={tuple(r[k] for k in ['guide','family','gene_a','gene_b','reference_gene']):r['triad_id'] for r in whole_links}
    associations=defaultdict(set)
    for r in domain_links:associations[r['triad_key']].add(lookup[tuple(r[k] for k in ['guide','family','gene_a','gene_b','reference_gene'])])
    assert len(associations)==26 and all(len(v)==1 for v in associations.values())
    dtriads={}
    for line in checked(dm,'triads.jsonl').open():
        r=json.loads(line)
        if r['triad_key'] in associations:dtriads[r['triad_key']]=r
    intervals={r[k] for r in dtriads.values() for k in ['interval_a','interval_b','interval_reference']}
    bounds={}
    for line in checked(base/'duplication-domain-inputs-20260926-v1','inputs.jsonl').open():
        r=json.loads(line)
        if r['interval_id'] in intervals and r['mask']=='full':bounds[r['interval_id']]=r
    wanted=set.union(*associations.values());models={(r[k+'_model'],int(r[k+'_version'])) for key,r in wtriads.items() if key in wanted for k in ['a','b','reference']}
    inputs={}
    for name in ['duplication-alignment-inputs-20260926-v1','duplication-reference-alignment-inputs-20260926-v1']:
        for line in checked(base/name,'inputs.jsonl').open():
            r=json.loads(line);key=(r['model_id'],r['version'])
            if key in models:inputs[key,r['mask']]=r
    @lru_cache(maxsize=80)
    def coordinate(key,mask):
        r=inputs[key,mask];seq,xyz,confidence=load_pdb(r)
        return {p:xyz[i] for i,p in enumerate(r['original_positions'])}
    maps={}
    with gzip.open(checked(wm,'common_residue_maps.jsonl.gz'),'rt') as f:
        for line in f:
            r=json.loads(line)
            if r['triad_id'] in wanted:maps[r['triad_id'],r['mask'],tuple(r['orders'])]=r
    assert len(maps)==13*16
    output=[];partition=[];maxerror=0.
    for dk in sorted(dtriads):
        wk=next(iter(associations[dk]));w=wtriads[wk];ds=dtriads[dk]
        modelkeys=[(w[k+'_model'],int(w[k+'_version'])) for k in ['a','b','reference']]
        bb=[bounds[ds['interval_'+k]] for k in ['a','b','reference']]
        assert [(b['model_id'],b['version']) for b in bb]==modelkeys
        for (tk,mask,orders),mapping in sorted(maps.items()):
            if tk!=wk:continue
            assert not any(mapping['edge_exclusions'])
            for definition,field in [('reference_common','reference_common_triples'),('cycle_consistent','cycle_consistent_triples')]:
                triples=mapping[field];inside=[];outside=[];mixed=[]
                for t in triples:
                    membership=[b['start']<=p<=b['end'] for b,p in zip(bb,t)]
                    (inside if all(membership) else outside if not any(membership) else mixed).append(t)
                # Independent set partition from interval memberships.
                bits=np.array([[b['start']<=t[i]<=b['end'] for i,b in enumerate(bb)] for t in triples],dtype=bool)
                assert [sum(bits.all(1)),sum(~bits.any(1)),sum(bits.any(1)&~bits.all(1))]==[len(inside),len(outside),len(mixed)]
                ident=dict(domain_triad=dk,whole_triad=wk,mask=mask,order_ab=orders[0],order_ar=orders[1],order_br=orders[2],mapping_definition=definition)
                partition.append(dict(ident,inside_triples=inside,outside_triples=outside,mixed_triples=mixed))
                coords=[coordinate(k,mask) for k in modelkeys]
                for pair,i,j in [('ab',0,1),('ar',0,2),('br',1,2)]:
                    r=dict(ident,pair=pair,inside_residues=len(inside),outside_residues=len(outside),mixed_residues=len(mixed),status='too_few_anchor_residues',anchor_rmsd='',outside_domain_anchored_rmsd='',outside_independent_fit_rmsd='',relative_anchor_curvature='')
                    if len(inside)>=3:
                        x,y=[np.array([coords[c][t[c]] for t in inside]) for c in [i,j]];g=geometry(x,y);r['relative_anchor_curvature']=g['relative_rotation_curvature']
                        if g['geometry_status']!='unique_at_numeric_tolerance':r['status']='degenerate_anchor'
                        else:
                            rot,trans=transform(x,y,'svd');qr,qt=transform(x,y,'quaternion');r['anchor_rmsd']=residual(x,y,rot,trans)
                            error=abs(r['anchor_rmsd']-residual(x,y,qr,qt));r['status']='no_outside_residues' if not outside else 'computed'
                            if outside:
                                ox,oy=[np.array([coords[c][t[c]] for t in outside]) for c in [i,j]]
                                r['outside_domain_anchored_rmsd']=residual(ox,oy,rot,trans)
                                error=max(error,abs(r['outside_domain_anchored_rmsd']-residual(ox,oy,qr,qt)))
                                if len(outside)>=3:
                                    rr,tt=transform(ox,oy,'svd');r['outside_independent_fit_rmsd']=residual(ox,oy,rr,tt)
                                    assert r['outside_independent_fit_rmsd']<=r['outside_domain_anchored_rmsd']+1e-8
                            assert error<1e-8;maxerror=max(maxerror,error)
                    output.append(r)
    assert len(partition)==832 and len(output)==2496
    out=base/'domain-anchored-displacement-20260927-v1';out.mkdir(exist_ok=False)
    path=out/'pair_displacements.tsv'
    with path.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(output[0]),delimiter='\t');writer.writeheader();writer.writerows(output)
    assert rows(path)==[{k:str(v) for k,v in r.items()} for r in output]
    mp=out/'residue_partitions.jsonl'
    with mp.open('w') as f:
        for r in partition:f.write(json.dumps(r,separators=(',',':'))+'\n')
    assert [json.loads(line) for line in mp.open()]==partition
    for path,digest in sources.items():assert sha(path)==digest
    receipt=dict(status='complete_domain_anchored_displacement_with_quaternion_checks',source_hashes=sources,script_sha256=sha(__file__),fixtures=fixtures,partitions=len(partition),pair_rows=len(output),statuses=dict(Counter(r['status'] for r in output)),maximum_quaternion_rmsd_difference=maxerror,artifacts={p.name:sha(p) for p in out.iterdir()},scope='Whole-protein residue correspondences partitioned by domain intervals: all three inside, all three outside, or mixed (retained but not fitted). Anchor fits require unique numerical geometry. Same domain transform applied to outside coordinates. Outside residuals include flexible regions and within/outside-domain changes; not proof of rigid interdomain rotation or biological significance. Native domain mappings are not substituted into whole-protein correspondence sets.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
