#!/usr/bin/env python3
"""Reconstruct every experimental quartet partition and verify all serialized fits."""
import csv,gzip,json,subprocess,time
from functools import lru_cache
from pathlib import Path
import numpy as np
import psutil
from screen_duplication_domain_alignment_coverage import sha

BASE=Path('results');ROOT=BASE/'experimental_structures/whole-domain-case-quartet-fits-20260927-v1'
OUTPUT=BASE/'experimental_structures/whole-domain-case-quartet-readback-20260927-v1'
ROLES=['a','b','reference'];PAIRS=[('ab',0,1),('ar',0,2),('br',1,2),('ae',0,3),('be',1,3),('re',2,3)]
ID=['triad_id','domain_triad','entity_id','experimental_model','label_asym_id','mask','context_indices_json']


def independent_fit(x,y):
    xm=x.mean(0);ym=y.mean(0);a=x-xm;b=y-ym;h=a.T@b
    xx,xy,xz=h[0];yx,yy,yz=h[1];zx,zy,zz=h[2]
    k=np.array([[xx+yy+zz,yz-zy,zx-xz,xy-yx],
                [yz-zy,xx-yy-zz,xy+yx,zx+xz],
                [zx-xz,xy+yx,-xx+yy-zz,yz+zy],
                [xy-yx,zx+xz,yz+zy,-xx-yy+zz]])
    values,vectors=np.linalg.eigh(k);w,i,j,l=vectors[:,-1]
    rotation=np.array([[1-2*(j*j+l*l),2*(i*j-w*l),2*(i*l+w*j)],
                       [2*(i*j+w*l),1-2*(i*i+l*l),2*(j*l-w*i)],
                       [2*(i*l-w*j),2*(j*l+w*i),1-2*(i*i+j*j)]]).T
    translation=ym-xm@rotation
    assert np.allclose(rotation.T@rotation,np.eye(3),atol=1e-12) and abs(np.linalg.det(rotation)-1)<1e-12
    u,s,vt=np.linalg.svd(h);curvature=s[1]+(1 if np.linalg.det(u@vt)>=0 else -1)*s[2]
    relative=curvature/s[0] if s[0]>0 else 0.
    unique=curvature>np.finfo(float).eps*max(len(x),3)*s[0]
    return rotation,translation,float(relative),bool(unique)


def rms(x,y,r,t):return float(np.sqrt(np.sum((x@r+t-y)**2)/len(x)))


def main():
    lp=Path('metadata/case_experimental_quartet_fits_launch_20260927.json');launch=json.loads(lp.read_text());lh=sha(lp)
    while True:
        try:
            p=psutil.Process(launch['pid'])
            if p.create_time()!=launch['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        print('Waiting for exact quartet fits',launch['pid'],flush=True);time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0') and sha(launch['cmdline'][1])==launch['script_sha256']
    rp=ROOT/'receipt.json';r=json.loads(rp.read_text());assert r['status']=='complete_case_experimental_quartet_geometric_fits'
    sources={**r['source_hashes'],str(lp):lh,str(rp):sha(rp)}
    for name,digest in r['artifacts'].items():sources[str(ROOT/name)]=digest
    for p,digest in sources.items():assert sha(p)==digest
    qr=BASE/'experimental_structures/whole-domain-case-sequence-quartets-20260927-v1'
    with gzip.open(qr/'quartet_residue_maps.jsonl.gz','rt') as f:quartets=list(map(json.loads,f))
    qlookup={(q['triad_id'],q['entity_id'],tuple(q['context_indices'])):q for q in quartets};assert len(qlookup)==len(quartets)
    structural=BASE/'structural_comparisons';case=structural/'whole-domain-case-dossiers-20260927-v1'
    def table(p):
        with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
    fields=['guide','family','gene_a','gene_b','reference_gene'];whole={tuple(x[k] for k in fields):x['triad_id'] for x in table(case/'whole_reference_links.tsv')}
    associations={}
    for x in table(case/'domain_reference_links.tsv'):
        value=whole[tuple(x[k] for k in fields)]
        if x['triad_key'] in associations:assert associations[x['triad_key']]==value
        associations[x['triad_key']]=value
    dt={}
    for line in (structural/'duplication-domain-common-residues-20260927-v1/triads.jsonl').open():
        x=json.loads(line)
        if x['triad_key'] in associations:dt[x['triad_key']]=x
    intervals={x['interval_'+role] for x in dt.values() for role in ROLES};bounds={}
    for line in (structural/'duplication-domain-inputs-20260926-v1/inputs.jsonl').open():
        x=json.loads(line)
        if x['mask']=='full' and x['interval_id'] in intervals:bounds[x['interval_id']]=x
    models={(m,v) for q in quartets for m,v in zip(q['model_ids'],q['model_versions'])};coords={}
    for name in ['duplication-alignment-inputs-20260926-v1','duplication-reference-alignment-inputs-20260926-v1']:
        for line in (structural/name/'inputs.jsonl').open():
            x=json.loads(line);key=(x['model_id'],x['version'])
            if key not in models:continue
            result={}
            if x['status']=='ready':
                assert sha(x['path'])==x['sha256']
                for line in Path(x['path']).read_text().splitlines():
                    if line.startswith('ATOM  '):
                        assert line[12:16].strip()=='CA';pos=int(line[22:26]);assert pos not in result
                        result[pos]=np.array([float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]])
                assert sorted(result)==x['original_positions']
            coords[(*key,x['mask'])]=result
    mapping=BASE/'experimental_structures/whole-domain-case-ca-mapping-20260927-v1';excluded=json.loads((mapping/'config.json').read_text())['deferred_entries']
    @lru_cache(maxsize=2)
    def experimental(entry):
        path=mapping/(entry+'.residues.tsv.gz');assert sha(path)==sources[str(path)];result={}
        with gzip.open(path,'rt') as f:
            for row in csv.DictReader(f,delimiter='\t'):
                key=(entry+'_'+row['entity_id'],row['model_number'],row['label_asym_id']);result.setdefault(key,{})
                if row['CA_status']=='unambiguous_full_occupancy_CA':
                    atom=json.loads(row['atom_records_json'])[0];result[key][int(row['label_seq_id'])]=np.array([float(atom[k]) for k in ['Cartn_x','Cartn_y','Cartn_z']])
        return result
    expected=set()
    for q in quartets:
        entry=q['entity_id'].rsplit('_',1)[0]
        if entry in excluded:continue
        for entity,model,chain in experimental(entry):
            if entity!=q['entity_id']:continue
            for domain,triad in associations.items():
                if triad!=q['triad_id']:continue
                for mask in ['full','plddt70']:expected.add((triad,domain,entity,model,chain,mask,tuple(q['context_indices'])))
    n_maps=n_fits=0;maximum=0.
    with (ROOT/'pair_fits.tsv').open() as f,gzip.open(ROOT/'observed_quartet_partitions.jsonl.gz','rt') as mf:
        fits=iter(csv.DictReader(f,delimiter='\t'))
        for m in map(json.loads,mf):
            indices=tuple(json.loads(m['context_indices_json']));key=tuple(m[k] for k in ID[:-1])+(indices,);assert key in expected;expected.remove(key);n_maps+=1
            q=qlookup[m['triad_id'],m['entity_id'],indices];modelkeys=list(zip(q['model_ids'],q['model_versions']));ec=experimental(m['entity_id'].rsplit('_',1)[0])[m['entity_id'],m['experimental_model'],m['label_asym_id']]
            cc=[coords[(*k,m['mask'])] for k in modelkeys]+[ec]
            common=[p for p in q['common_positions'] if all(p[i] in cc[i] for i in range(4))];assert m['common_positions']==common
            bb=[bounds[dt[m['domain_triad']]['interval_'+role]] for role in ROLES];assert [(b['model_id'],b['version']) for b in bb]==modelkeys
            memberships=[sum(b['start']<=p[i]<=b['end'] for i,b in enumerate(bb)) for p in common]
            inside=[p for p,c in zip(common,memberships) if c==3];outside=[p for p,c in zip(common,memberships) if c==0];mixed=[p for p,c in zip(common,memberships) if c in [1,2]]
            assert [m[k] for k in ['inside_positions','outside_positions','mixed_positions']]==[inside,outside,mixed]
            assert m['original_query_lengths']==q['query_lengths'] and m['domain_lengths']==[b['end']-b['start']+1 for b in bb] and m['outside_lengths']==[b['original_length']-(b['end']-b['start']+1) for b in bb]
            for region,positions in [('whole',common),('domain',inside),('outside',outside)]:
                for pair,i,j in PAIRS:
                    row=next(fits);n_fits+=1;assert all(row[k]==m[k] for k in ID) and row['region']==region and row['pair']==pair and int(row['residues'])==len(positions)
                    if len(positions)<3:
                        assert row['status']=='too_few_residues' and all(row[k]=='' for k in ['rmsd','relative_rotation_curvature','outside_under_domain_transform_rmsd']);continue
                    x,y=[np.array([cc[k][p[k]] for p in positions]) for k in [i,j]];rot,trans,curvature,unique=independent_fit(x,y)
                    assert abs(float(row['relative_rotation_curvature'])-curvature)<1e-10
                    if not unique:assert row['status']=='degenerate_rotation' and row['rmsd']=='' and row['outside_under_domain_transform_rmsd']=='';continue
                    assert row['status']=='computed';error=abs(float(row['rmsd'])-rms(x,y,rot,trans))
                    if region=='domain' and outside:
                        ox,oy=[np.array([cc[k][p[k]] for p in outside]) for k in [i,j]];error=max(error,abs(float(row['outside_under_domain_transform_rmsd'])-rms(ox,oy,rot,trans)))
                    else:assert row['outside_under_domain_transform_rmsd']==''
                    assert error<1e-8;maximum=max(maximum,error)
            if n_maps%100==0:print(n_maps,'partitions verified',flush=True)
        assert next(fits,None) is None
    assert not expected and n_maps==r['partition_maps'] and n_fits==r['pair_fit_rows']==18*n_maps
    for p,digest in sources.items():assert sha(p)==digest
    OUTPUT.mkdir(exist_ok=False);result=dict(status='complete_full_experimental_quartet_readback',terminal_state=state,source_hashes=sources,script_sha256=sha(__file__),partition_maps=n_maps,pair_fit_rows=n_fits,maximum_rmsd_disagreement=maximum,scope='Full expected combination set, residue availability, domain partition and denominator reconstruction; every serialized proper-rotation fit and outside-under-domain-transform RMSD checked using separately implemented quaternion matrix formulas. Shared numeric linear algebra library and source annotations remain dependencies. No coverage qualification or biological conclusion.')
    (OUTPUT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
