#!/usr/bin/env python3
"""Describe both directed PAE blocks between mapped domain anchors and outside residues."""
import csv,gzip,json,math
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha


def main():
    base=Path('results/structural_comparisons');sources={}
    auditpath=Path('metadata/whole_domain_case_pae_readback_20260927.json');audit=json.loads(auditpath.read_text())
    assert audit['status']=='passed_full_case_pae_retrieval_readback' and audit['models_verified']==39 and audit['models_failed']==0
    sources[str(auditpath)]=sha(auditpath)
    for path,digest in audit['source_hashes'].items():assert sha(path)==digest
    def checked(root,name):
        rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/name;assert sha(p)==r['artifacts'][name]
        sources[str(rp)]=sha(rp);sources[str(p)]=sha(p);return p
    manifest=json.loads(checked(base/'whole-domain-case-pae-20260927-v1','pae_manifest.json').read_text())
    matrices={}
    for r in manifest:
        p=Path(r['path']);assert sha(p)==r['gzip_sha256'];sources[str(p)]=sha(p)
        matrices[r['model_id'],r['version']]=np.asarray(json.loads(gzip.decompress(p.read_bytes()))[0]['predicted_aligned_error'],dtype=float)
    with checked(base/'whole-protein-common-residues-20260927-v1','model_triads.tsv').open() as f:models={r['triad_id']:r for r in csv.DictReader(f,delimiter='\t')}
    partitions=[json.loads(line) for line in checked(base/'domain-anchored-displacement-20260927-v1','residue_partitions.jsonl').open()]
    @lru_cache(maxsize=None)
    def statistics(model,version,ii,jj,diagonal):
        matrix=matrices[model,version];i=np.array(ii)-1;j=np.array(jj)-1
        assert len(set(ii))==len(ii) and len(set(jj))==len(jj) and min(ii+jj)>=1 and max(ii+jj)<=len(matrix)
        block=matrix[np.ix_(i,j)]
        values=block[~np.eye(len(i),dtype=bool)] if diagonal else block.ravel()
        # Independent scalar indexing, sum, rank interpolation and threshold counts.
        scalar=[float(matrix[a-1,b-1]) for a in ii for b in jj if not diagonal or a!=b]
        ordered=sorted(scalar);assert len(values)==len(scalar)>0
        r=dict(entries=len(values),mean=float(np.mean(values)),median=float(np.quantile(values,.5)),p90=float(np.quantile(values,.9)),maximum=float(np.max(values)))
        assert math.isclose(r['mean'],math.fsum(scalar)/len(scalar),rel_tol=1e-12,abs_tol=1e-12) and r['maximum']==ordered[-1]
        for label,q in [('median',.5),('p90',.9)]:
            position=(len(ordered)-1)*q;lo=int(position);hi=min(lo+1,len(ordered)-1)
            assert math.isclose(r[label],ordered[lo]+(position-lo)*(ordered[hi]-ordered[lo]),rel_tol=1e-12,abs_tol=1e-12)
        for threshold in [5,10,15]:
            r[f'fraction_le_{threshold}']=float(np.mean(values<=threshold))
            assert r[f'fraction_le_{threshold}']==sum(v<=threshold for v in scalar)/len(scalar)
        return r
    output=[]
    for p in partitions:
        t=models[p['whole_triad']]
        for column,role in enumerate(['a','b','reference']):
            inside=tuple(sorted(x[column] for x in p['inside_triples']));outside=tuple(sorted(x[column] for x in p['outside_triples']))
            assert not set(inside)&set(outside)
            for label,ii,jj,diagonal in [('inside_inside_offdiagonal',inside,inside,True),('row_inside_col_outside',inside,outside,False),('row_outside_col_inside',outside,inside,False)]:
                ident={k:p[k] for k in ['domain_triad','whole_triad','mask','order_ab','order_ar','order_br','mapping_definition']}
                output.append(dict(ident,role=role,model_id=t[role+'_model'],version=int(t[role+'_version']),block=label,**statistics(t[role+'_model'],int(t[role+'_version']),ii,jj,diagonal)))
    assert len(output)==832*3*3
    links=pd.read_csv(checked(base/'whole-domain-case-dossiers-20260927-v1','domain_reference_links.tsv'),sep='\t')
    case=['family','gene_a','gene_b','pfam_accession'];relation=links[case+['triad_key']].drop_duplicates().rename(columns={'triad_key':'domain_triad'})
    data=pd.DataFrame(output);linked=relation.merge(data,on='domain_triad',validate='one_to_many');assert len(linked)==len(data)
    summary=[]
    for key,group in linked.groupby(case+['role','block'],sort=True):
        r=dict(zip(case+['role','block'],key),alternatives=len(group))
        for field in ['mean','median','p90','maximum','fraction_le_5','fraction_le_10','fraction_le_15']:
            r[field+'_min']=float(group[field].min());r[field+'_max']=float(group[field].max())
        summary.append(r)
    assert len(summary)==117
    out=base/'case-regional-pae-20260927-v1';out.mkdir(exist_ok=False);artifacts={}
    for name,frame in [('all_regional_pae.tsv',linked),('case_regional_ranges.tsv',pd.DataFrame(summary))]:
        p=out/name;frame.to_csv(p,sep='\t',index=False);back=pd.read_csv(p,sep='\t')
        pd.testing.assert_frame_equal(back,frame.reset_index(drop=True),check_dtype=False,rtol=1e-12,atol=1e-12);artifacts[name]=sha(p)
    for path,digest in sources.items():assert sha(path)==digest
    receipt=dict(status='complete_case_regional_pae_with_scalar_readback',source_hashes=sources,script_sha256=sha(__file__),regional_rows=len(output),case_role_block_rows=len(summary),unique_region_statistics=statistics.cache_info().currsize,artifacts=artifacts,
        scope='Raw matrix row/column directions retained separately; domain internal diagonal removed. All mapped anchor/outside residue sets, masks, boundaries and orders retained. Full scalar indexing and quantile interpolation verify each unique region. 5/10/15-A fractions are descriptive sensitivity thresholds, not calibrated acceptance criteria or structural-contrast confidence intervals.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
