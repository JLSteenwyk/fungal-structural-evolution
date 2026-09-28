#!/usr/bin/env python3
"""Compare conditional marginals at exact projected coordinate matches."""
from collections import defaultdict
import csv,gzip,hashlib,json
from pathlib import Path
import numpy as np

ROOT=Path('results/ancestral')


def main():
    pins={}
    def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    def read(path):
        path=Path(path);pins[str(path)]=sha(path);return json.loads(path.read_text())
    corr=ROOT/'whole-domain-alignment-correspondence-20260927-v1'
    cr=read(corr/'receipt.json')
    for name,h in cr['artifacts'].items():assert sha(corr/name)==h
    with (corr/'comparison_summary.tsv').open() as h:summaries=list(csv.DictReader(h,delimiter='\t'))
    coordinates=defaultdict(list)
    with gzip.open(corr/'column_correspondence.tsv.gz','rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            if r['status']=='matched_projected_column':coordinates[tuple(r[k] for k in ['family','boundary','whole_method','domain_method'])].append(r)
    mp=ROOT/'case-local-trees-20260927-v1/ancestral_node_mapping.tsv';pins[str(mp)]=sha(mp)
    with mp.open() as h:mappings={(r['family'],r['dataset'],int(r['level'])):r for r in csv.DictReader(h,delimiter='\t') if r['guide']=='profile'}
    folders={'whole':ROOT/'refined-whole-ancestors-20260927-v2','domain':ROOT/'refined-domain-ancestors-20260927-v1'}
    receipts={label:read(folder/'receipt.json') for label,folder in folders.items()}
    for label in folders:
        completion=read('metadata/refined_'+label+'_posterior_audit_completed_20260927.json')
        path=Path(completion['completed_receipt_path'])
        assert sha(path)==completion['completed_receipt_sha256'];read(path)
    cache={}
    def posterior(label,jid):
        key=(label,jid)
        if key not in cache:
            path=folders[label]/(jid+'.npz');assert sha(path)==receipts[label]['artifacts'][path.name]
            with np.load(path) as z:
                assert z['levels'].tolist()==[0,1,2]
                assert ''.join(z['amino_acids'].tolist())=='ARNDCQEGHILKMFPSTWYV'
                cache[key]=z['posterior'].copy()
        return cache[key]
    out=ROOT/'refined-whole-domain-probability-comparison-20260927-v1';out.mkdir(exist_ok=False)
    fields=['family','boundary','whole_method','domain_method','model','alpha_min','level','source_node','whole_job','domain_job','whole_column','domain_column','coordinate_signature_sha256','whole_map','domain_map','whole_map_probability','domain_map_probability','total_variation','opposing_at_least_090']
    rows=[];aa='ARNDCQEGHILKMFPSTWYV';total=0;changed=0;opposed=0;maximum=0.
    with gzip.open(out/'matched_site_comparisons.tsv.gz','wt') as handle:
        writer=csv.DictWriter(handle,fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for summary in summaries:
            fam=summary['family'];boundary=summary['boundary'];wm=summary['whole_method'];dm=summary['domain_method']
            cols=coordinates[fam,boundary,wm,dm];assert len(cols)==int(summary['exact_shared_columns'])
            wi=np.array([int(r['whole_column'])-1 for r in cols],dtype=int);di=np.array([int(r['domain_column'])-1 for r in cols],dtype=int)
            retained=set(json.loads(mappings[fam,'domain',3]['retained_set_json']))
            for model in ['LG','WAG','JTT']:
                for bound in ['0.02','0.005']:
                    wj=f'{fam}-whole-{wm}-{model}-refine-amin{bound}';dj=f'{fam}-{boundary}-{dm}-{model}-refine-amin{bound}'
                    wp=posterior('whole',wj);dp=posterior('domain',dj)
                    for level in range(3):
                        wn=mappings[fam,'whole',level];dn=mappings[fam,'domain',level]
                        assert wn['source_node']==dn['source_node'] and wn['unique_output_node']==dn['unique_output_node']=='1'
                        ws=set(json.loads(wn['retained_set_json']));ds=set(json.loads(dn['retained_set_json']));assert ws&retained==ds
                        x=wp[level,wi];y=dp[level,di];tv=np.abs(x-y).sum(axis=1)/2
                        xm=x.argmax(axis=1);ym=y.argmax(axis=1);xp=x.max(axis=1);yp=y.max(axis=1)
                        different=xm!=ym;high=different&(xp>=.9)&(yp>=.9)
                        common=dict(family=fam,boundary=boundary,whole_method=wm,domain_method=dm,model=model,alpha_min=bound,level=level,source_node=wn['source_node'],whole_job=wj,domain_job=dj)
                        for i,c in enumerate(cols):
                            writer.writerow(dict(common,whole_column=c['whole_column'],domain_column=c['domain_column'],coordinate_signature_sha256=c['coordinate_signature_sha256'],whole_map=aa[xm[i]],domain_map=aa[ym[i]],whole_map_probability=xp[i],domain_map_probability=yp[i],total_variation=tv[i],opposing_at_least_090=bool(high[i])))
                        rows.append(dict(common,whole_descendants=len(ws),domain_descendants=len(ds),matched_sites=len(cols),only_whole_projection=int(summary['only_whole_projection']),only_domain_alignment=int(summary['only_domain_alignment']),map_disagreements=int(different.sum()),opposing_at_least_090=int(high.sum()),maximum_total_variation=float(tv.max()) if len(tv) else None,total_variation_sum=float(tv.sum())))
                        total+=len(cols);changed+=int(different.sum());opposed+=int(high.sum());maximum=max(maximum,float(tv.max()) if len(tv) else 0.)
    assert len(rows)==1872 and len(cache)==468
    with (out/'node_summary.tsv').open('w') as h:
        w=csv.DictWriter(h,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    result=dict(status='complete_refined_whole_domain_conditional_comparison',alignment_comparisons=104,fit_pairs=624,node_comparisons=len(rows),matched_node_sites=total,map_disagreements=changed,opposing_at_least_090=opposed,maximum_total_variation=maximum,pins=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Exact shared coordinate signatures only; unmatched counts retained per node and full coordinate union retained in source. Source-node identity and restriction of descendants checked. Differences combine protein membership, flanking context, alignment and parameter estimation; they do not isolate a causal context effect. Dependent comparisons are not evolutionary events or posterior mixture weights. Conditional on residue existence and assumed local tree; indel uncertainty unresolved.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['pins','artifacts']}))


if __name__=='__main__':main()
