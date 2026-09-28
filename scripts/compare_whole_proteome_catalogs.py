"""Compare two validated catalogs without equating download growth with coverage."""
import argparse
import json
from pathlib import Path
import pandas as pd
from readback_whole_proteome_catalog import sha

KEYS=['taxon_id','protein_id']
MODEL=['model_id','version','model_path']


def compare_links(old,new):
    for frame in [old,new]:
        assert not frame.duplicated(KEYS).any()
        assert frame[KEYS+MODEL+['sequence_sha256']].notna().all().all()
    joined=old.merge(new,on=KEYS,how='outer',suffixes=('_old','_new'),indicator=True,validate='one_to_one')
    common=joined['_merge']=='both'
    assert (joined.loc[common,'sequence_sha256_old']==joined.loc[common,'sequence_sha256_new']).all(), 'Protein sequence universe changed'
    joined['disposition']='unchanged_model'
    joined.loc[joined['_merge']=='left_only','disposition']='lost_catalog_link'
    joined.loc[joined['_merge']=='right_only','disposition']='new_catalog_link'
    changed=pd.Series(False,index=joined.index)
    for col in MODEL:changed |= joined[col+'_old']!=joined[col+'_new']
    joined.loc[common & changed,'disposition']='changed_selected_model'
    return joined.drop(columns='_merge')


def load_catalog(root,audit):
    r=json.loads((root/'receipt.json').read_text());a=json.loads(audit.read_text())
    assert r['status']=='complete_whole_representative_proteome_exact_sequence_catalog'
    assert a['status']=='passed_full_proteome_sequence_and_model_selection_readback'
    assert a['producer_receipt_sha256']==sha(root/'receipt.json')
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    links=pd.read_csv(root/'protein_model_links.tsv',sep='\t',dtype=str,keep_default_na=False)
    coverage=pd.read_csv(root/'taxon_coverage.tsv',sep='\t')
    assert len(links)==r['proteins_linked'] and len(coverage)==r['taxa']
    return r,links,coverage


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    p=json.loads(args.plan.read_text())
    for name,h in p['pins'].items():assert sha(name)==h
    roots=[Path(p[k]) for k in ['old','new']]
    old,ol,oc=load_catalog(roots[0],Path(p['old_audit']))
    new,nl,nc=load_catalog(roots[1],Path(p['new_audit']))
    assert old['proteins_screened']==new['proteins_screened']==5815847
    assert old['taxa']==new['taxa']==526
    for frame in [oc,nc]:assert not frame.taxon_id.duplicated().any()
    stable=['taxon_id','species_name','study_role','representative_proteins']
    pd.testing.assert_frame_equal(oc[stable].sort_values('taxon_id').reset_index(drop=True),nc[stable].sort_values('taxon_id').reset_index(drop=True))
    joined=compare_links(ol,nl)
    counts=pd.crosstab(joined.taxon_id,joined.disposition).reindex(columns=['new_catalog_link','lost_catalog_link','changed_selected_model','unchanged_model'],fill_value=0)
    table=oc.merge(nc,on=stable,validate='one_to_one',suffixes=('_old','_new')).merge(counts,on='taxon_id',how='left',validate='one_to_one')
    for col in counts.columns:table[col]=table[col].fillna(0).astype(int)
    assert (table.proteins_with_model_old==table.lost_catalog_link+table.changed_selected_model+table.unchanged_model).all()
    assert (table.proteins_with_model_new==table.new_catalog_link+table.changed_selected_model+table.unchanged_model).all()
    table['unlinked_in_both']=table.representative_proteins-table[counts.columns].sum(axis=1)
    assert (table.unlinked_in_both>=0).all()
    for when in ['old','new']:
        table['coverage_fraction_'+when]=table['proteins_with_model_'+when]/table.representative_proteins
    table['coverage_fraction_change']=table.coverage_fraction_new-table.coverage_fraction_old
    out=Path(p['output']);out.mkdir(parents=True,exist_ok=False)
    joined.sort_values(KEYS).to_csv(out/'protein_link_dispositions.tsv',sep='\t',index=False)
    table.sort_values('taxon_id').to_csv(out/'taxon_coverage_change.tsv',sep='\t',index=False)
    receipt=dict(status='complete_validated_catalog_comparison',plan_sha256=sha(args.plan),
        old_links=len(ol),new_links=len(nl),net_link_change=len(nl)-len(ol),
        dispositions=joined.disposition.value_counts().to_dict(),unlinked_in_both=int(table.unlinked_in_both.sum()),
        sources={str(path):sha(path) for path in [roots[0]/'receipt.json',roots[1]/'receipt.json',Path(p['old_audit']),Path(p['new_audit'])]},
        artifacts={f.name:sha(f) for f in out.iterdir()},
        scope='Exact-sequence catalog coverage change in a fixed representative-protein universe. Lost links and model replacements retained. '
              'Newly covered proteins are not necessarily newly inferred structures; no new confidence qualification or evolutionary inference.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))

if __name__=='__main__':main()
