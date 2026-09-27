#!/usr/bin/env python3
"""Independently check every robustness summary through dataframe grouping of source fits."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();p=json.loads(a.plan.read_text());root=Path(p['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
    assert r['status']=='complete_common_core_robustness_pending_readback' and r['plan_sha256']==sha(a.plan)
    for path,h in p['pins'].items():assert sha(path)==h,path
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    assert sha(Path(p['fits'])/'receipt.json')==r['fit_receipt_sha256'] and sha(p['fit_readback'])==r['fit_readback_sha256']
    fields=['triad_key','mask','mapping_definition','order_ab','order_ar','order_br','rmsd_ar_minus_br']+[s+'_pass' for s in p['screens']]
    fits=pd.read_csv(Path(p['fits'])/'common_residue_fits.tsv',sep='\t',usecols=fields)
    assert not fits.duplicated(fields[:6]).any()
    assert set(fits['mask'])=={'full','plddt70'} and set(fits['mapping_definition'])=={'reference_common','cycle_consistent'}
    for col in ['order_ab','order_ar','order_br']:assert set(fits[col])=={0,1}
    assert fits.groupby('triad_key').size().eq(32).all()
    actual=pd.read_csv(root/'robustness.tsv',sep='\t',keep_default_na=True)
    keys=['triad_key','scope','mask','mapping_definition','screen'];assert not actual.duplicated(keys).any()
    blocks=[]
    for screen in p['screens']:
        assert set(fits[screen+'_pass'])<={0,1}
        eligible=fits.loc[fits[screen+'_pass'].eq(1)]
        assert np.isfinite(eligible['rmsd_ar_minus_br']).all()
        for scope,group,n in [('input_orders',['triad_key','mask','mapping_definition'],8),('all_alternatives',['triad_key'],32)]:
            universe=fits.groupby(group).size();assert universe.eq(n).all()
            values=eligible.groupby(group)['rmsd_ar_minus_br'].agg(['size','min','max']).reindex(universe.index)
            values['eligible_alternatives']=values['size'].fillna(0).astype(int);values['expected_alternatives']=n;values['complete']=values['eligible_alternatives'].eq(n).astype(int)
            values['eligible_contrast_min']=values['min'];values['eligible_contrast_max']=values['max'];values['eligible_contrast_range']=values['max']-values['min']
            for margin in p['margins']:
                m=margin['angstrom'];complete=values['complete'].eq(1)
                values[margin['id']]=np.select([~complete,values['min'].gt(m),values['max'].lt(-m),values['min'].lt(-m)&values['max'].gt(m)],['incomplete','a_farther_from_reference','b_farther_from_reference','changes_direction_across_alternatives'],default='touches_or_within_margin_band')
            values=values.drop(columns=['size','min','max']).reset_index();values['scope']=scope;values['screen']=screen
            if scope=='all_alternatives':values['mask']='all';values['mapping_definition']='all'
            blocks.append(values)
    expected=pd.concat(blocks,ignore_index=True).set_index(keys).sort_index();actual=actual.set_index(keys).sort_index()
    assert expected.index.equals(actual.index)
    for col in expected:
        if col.startswith('eligible_contrast_'):assert np.allclose(expected[col],actual[col],rtol=1e-12,atol=1e-12,equal_nan=True),col
        else:assert expected[col].equals(actual[col]),col
    counts={}
    for margin in p['margins']:
        groups=actual.reset_index().groupby(['scope','mask','mapping_definition','screen',margin['id']]).size()
        for k,n in groups.items():counts['|'.join([*k[:4],margin['id'],k[4]])]=int(n)
    assert counts==r['classification_counts'] and len(actual)==r['summary_rows'] and len(fits)==r['source_fit_rows']
    selected=actual.reset_index().query("scope == 'all_alternatives' and screen == 'n30_c70'")
    selected_counts={m['id']:selected[m['id']].value_counts().to_dict() for m in p['margins']}
    result=dict(status='passed_full_common_core_robustness_readback',producer_receipt_sha256=sha(rp),plan_sha256=sha(a.plan),checker_sha256=sha(__file__),summary_rows_checked=len(actual),source_rows_checked=len(fits),triads=int(fits['triad_key'].nunique()),n30_c70_all_alternatives_classifications=selected_counts,scope='Every source alternative grid and all summary keys, counts, extrema, ranges, completeness flags and direction classifications independently recomputed using dataframe group operations. Margins descriptive only; triads are not independent events and no biological significance is inferred.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
