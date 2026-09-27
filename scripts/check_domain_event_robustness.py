#!/usr/bin/env python3
"""Independently reconstruct event/domain robustness by relational merges and group operations."""
import argparse,hashlib,json,sqlite3
from pathlib import Path
import numpy as np
import pandas as pd

E=['guide','family','gene_node','gene_a','gene_b']
S=['taxon_id','chosen_reference_gene','duplicate_pair_sequence_distance','sequence_signed_difference','sequence_normalized_contrast','sequence_direction','sequence_covariate_status']

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def check_table(path,expected,keys,numeric=()):
    actual=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False)
    assert not actual.duplicated(keys).any() and not expected.duplicated(keys).any()
    actual=actual.set_index(keys).sort_index();expected=expected.set_index(keys).sort_index()
    assert actual.index.equals(expected.index) and set(actual.columns)==set(expected.columns)
    for c in expected:
        if c in numeric:
            values=pd.to_numeric(actual[c].replace('',np.nan));assert np.allclose(values,expected[c],rtol=1e-12,atol=1e-12,equal_nan=True),c
        else:assert actual[c].equals(expected[c].astype(str)),c
    return actual.reset_index()

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();p=json.loads(a.plan.read_text());root=Path(p['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
    assert r['status']=='complete_domain_event_robustness_pending_readback' and r['plan_sha256']==sha(a.plan)
    for path,h in p['pins'].items():assert sha(path)==h,path
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    db=sqlite3.connect(f"file:{p['integration']}/domain_triads.sqlite?mode=ro",uri=True)
    triads=pd.read_sql_query('SELECT '+','.join(E+['reference_gene','policy'])+' FROM triads',db)
    sequence=pd.read_sql_query('SELECT '+','.join(E+S)+' FROM sequence_covariates',db);db.close()
    assert not triads.duplicated(E+['reference_gene','policy']).any() and not sequence.duplicated(E).any()
    universe=triads.groupby(E).agg(reference_choices=('reference_gene','nunique'),expected_reference_policy_boundary_units=('policy','size')).reset_index()
    universe['expected_reference_policy_boundary_units']*=2
    assert universe['expected_reference_policy_boundary_units'].eq(universe['reference_choices']*8).all()
    links=pd.read_csv(Path(p['mappings'])/'event_domain_links.tsv',sep='\t',dtype=str,keep_default_na=False)
    units=E+['reference_gene','policy','boundary','pfam_accession'];assert not links.duplicated(units).any()
    source=pd.read_csv(Path(p['robustness'])/'robustness.tsv',sep='\t')
    source=source.loc[source['scope'].eq('all_alternatives'),['triad_key','screen','complete','eligible_contrast_min','eligible_contrast_max']]
    assert not source.duplicated(['triad_key','screen']).any()
    linked=links.merge(source,on='triad_key',how='left',validate='many_to_many')
    assert len(linked)==len(links)*len(p['screens']) and not linked.duplicated(units+['screen']).any() and linked['screen'].notna().all()
    keys=E+['pfam_accession','screen']
    expected=linked.groupby(keys).agg(observed_units=('triad_key','size'),complete_units=('complete','sum'),eligible_contrast_min=('eligible_contrast_min','min'),eligible_contrast_max=('eligible_contrast_max','max')).reset_index()
    expected=expected.merge(universe,on=E,validate='many_to_one').merge(sequence,on=E,validate='many_to_one')
    expected['missing_units']=expected['expected_reference_policy_boundary_units']-expected['observed_units'];assert expected['missing_units'].ge(0).all()
    expected['complete']=expected['complete_units'].eq(expected['expected_reference_policy_boundary_units']).astype(int)
    for margin in p['margins']:
        m=margin['angstrom'];lo=expected['eligible_contrast_min'];hi=expected['eligible_contrast_max']
        expected[margin['id']]=np.select([expected['complete'].eq(0),lo.gt(m),hi.lt(-m),lo.lt(-m)&hi.gt(m)],['incomplete','a_farther_from_reference','b_farther_from_reference','changes_direction_across_alternatives'],default='touches_or_within_margin_band')
    actual=check_table(root/'event_domain_robustness.tsv',expected,keys,['eligible_contrast_min','eligible_contrast_max'])
    domains=links.groupby(E)['pfam_accession'].nunique().rename('domain_candidates').reset_index()
    availability=universe[E+['reference_choices']].merge(domains,on=E,how='left',validate='one_to_one').merge(sequence,on=E,validate='one_to_one')
    availability['domain_candidates']=availability['domain_candidates'].fillna(0).astype(int)
    availability['availability']=np.where(availability['domain_candidates'].gt(0),'domain_candidates_present','no_common_domain_candidate')
    check_table(root/'event_availability.tsv',availability,E)
    assert len(availability)==r['event_universe'] and len(expected)==r['event_domain_screen_rows'] and int(availability['domain_candidates'].sum())==r['event_domain_combinations']
    assert int(availability['domain_candidates'].gt(0).sum())==r['events_with_candidates']
    counts={}
    for margin in p['margins']:
        for (guide,screen,label),n in expected.groupby(['guide','screen',margin['id']]).size().items():counts['|'.join([guide,screen,margin['id'],label])]=int(n)
    assert counts==r['classification_counts']
    selected=[]
    for margin in p['margins']:
        for (guide,label),n in expected.loc[expected['screen'].eq('n30_c70')].groupby(['guide',margin['id']]).size().items():selected.append(dict(guide=guide,margin=margin['id'],classification=label,event_domain_combinations=int(n)))
    result=dict(status='passed_full_domain_event_robustness_readback',producer_receipt_sha256=sha(rp),plan_sha256=sha(a.plan),checker_sha256=sha(__file__),events_checked=len(availability),event_domain_screen_rows_checked=len(expected),event_domain_combinations=r['event_domain_combinations'],events_with_domain_candidates=r['events_with_candidates'],n30_c70_classifications=selected,scope='All event availability rows, original sequence covariates, complete tied-reference/policy/boundary grids, interval-summary joins, extrema, completeness and direction labels independently reconstructed. Event-domain rows and guides remain dependent; no biological significance or ancestral-direction inference.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
