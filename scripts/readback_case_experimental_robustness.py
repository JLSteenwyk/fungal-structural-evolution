#!/usr/bin/env python3
"""Recompute all case robustness ranges, coverage decisions and summary counts."""
import csv,json
from pathlib import Path
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha

BASE=Path('results/experimental_structures')
ROOT=BASE/'whole-domain-case-reference-robustness-20260927-v1'
SOURCE=BASE/'whole-domain-case-quartet-contrasts-20260927-v1'
UNIT=['triad_id','entity_id','experimental_model','label_asym_id','context_indices_json']
CASE=['family','gene_a','gene_b','pfam_accession']


def main():
    rp=ROOT/'receipt.json';r=json.loads(rp.read_text());assert r['status']=='complete_descriptive_case_reference_robustness'
    for p,digest in r['source_hashes'].items():assert sha(p)==digest
    for name,digest in r['artifacts'].items():assert sha(ROOT/name)==digest
    load=lambda p:pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)
    contrasts=load(SOURCE/'paired_reference_contrasts.tsv');coverage=load(SOURCE/'coverage_screens.tsv');data=load(ROOT/'reference_unit_sensitivity.tsv');summary=load(ROOT/'case_summary.tsv')
    keys=UNIT+['metric'];numerics=['fungal_reference_ar_minus_br','experimental_reference_ae_minus_be']
    for k in numerics:contrasts[k]=pd.to_numeric(contrasts[k],errors='coerce')
    grouped=contrasts.groupby(keys,dropna=False);counts=grouped.size();assert (counts==4).all()
    expected=grouped[numerics].agg(['min','max','count'])
    exact={}
    for key,row in expected.iterrows():exact[key]=row
    idcols=UNIT+['domain_triad','mask']
    projected=contrasts[idcols+['metric','coverage_requirement']].merge(coverage,on=idcols,validate='many_to_many')
    assert len(projected)==len(contrasts)*6
    projected['pass']=projected.apply(lambda row:row[row['coverage_requirement']+'_pass']=='1',axis=1)
    allpass=projected.groupby(keys+['screen'])['pass'].all().to_dict()
    seen=set()
    for record in data.to_dict('records'):
        key=tuple(record[k] for k in keys);screen=record['screen'];margin=float(record['margin_angstrom']);identity=key+(screen,margin)
        assert identity not in seen;seen.add(identity);original=exact[key]
        numeric=all(original[k,'count']==4 for k in numerics);coverage_pass=bool(allpass[key+(screen,)])
        assert int(record['required_variants'])==4 and int(record['all_numeric_pass'])==numeric and int(record['all_coverage_pass'])==coverage_pass
        qualified=numeric and coverage_pass;assert int(record['qualified'])==qualified
        directions=[]
        for prefix,col in zip(['fungal','experimental'],numerics):
            for stat in ['min','max']:
                value=original[col,stat];reported=record[prefix+'_'+stat]
                if pd.isna(value):assert reported==''
                else:assert abs(float(reported)-value)<1e-12
            status='not_qualified'
            if qualified:
                low,high=original[col,'min'],original[col,'max']
                status='positive' if low>margin+1e-8 else 'negative' if high<-margin-1e-8 else 'variable_or_within_margin'
            assert record[prefix+'_direction']==status;directions.append(status)
        status='not_qualified' if not qualified else 'same_direction' if directions[0]==directions[1] and directions[0] in ['positive','negative'] else 'opposite_direction' if set(directions)=={'positive','negative'} else 'variable_or_within_margin'
        assert record['reference_agreement']==status
    expected_keys={key+(f'n{n}_c{c}',margin) for key in exact for n in [30,50] for c in [50,70,90] for margin in [0.,.01,.1]};assert seen==expected_keys
    groupcols=CASE+['metric','screen','margin_angstrom'];lookup={key:frame for key,frame in data.groupby(groupcols)};seen_summaries=set()
    for row in summary.to_dict('records'):
        key=tuple(row[k] for k in groupcols);assert key not in seen_summaries;seen_summaries.add(key)
        frame=lookup.get(key,data.iloc[:0]);eligible=frame[frame.qualified=='1']
        assert int(row['candidate_units'])==len(frame) and int(row['qualified_units'])==len(eligible)
        for field,column in [('qualified_entities','entity_id'),('qualified_sequences','sequence_sha256'),('qualified_dependency_components','dependency_component')]:assert int(row[field])==eligible[column].nunique()
        assert int(row['qualified_entries'])==eligible.entity_id.str.rsplit('_',n=1).str[0].nunique()
        for status in ['same_direction','opposite_direction','variable_or_within_margin']:
            chosen=eligible[eligible.reference_agreement==status];assert int(row[status+'_units'])==len(chosen) and int(row[status+'_entities'])==chosen.entity_id.nunique()
    assert len(summary)==936 and len(data)==r['unit_sensitivity_rows'] and len(summary)==r['case_summary_rows']
    cases=load(Path('results/structural_comparisons/whole-domain-case-dossiers-20260927-v1/case_dossiers.tsv'))
    desired={tuple(row[k] for k in CASE)+(metric,f'n{n}_c{c}',str(margin)) for row in cases.to_dict('records') for metric in ['whole','domain','outside_independent','outside_domain_anchored'] for n in [30,50] for c in [50,70,90] for margin in [0.,.01,.1]}
    assert seen_summaries==desired
    result=dict(status='complete_full_case_reference_robustness_readback',source_receipt_sha256=sha(rp),checker_sha256=sha(__file__),source_contrast_receipt_sha256=sha(SOURCE/'receipt.json'),unit_rows=len(data),case_summary_rows=len(summary),scope='Every range/numeric count recomputed with pandas grouping from source contrasts; every coverage decision independently joined from source screens. All margins/directions and complete unit grid verified. All zero-inclusive case summaries recomputed; no inferential validation implied.')
    out=BASE/'whole-domain-case-reference-robustness-readback-20260927-v1';out.mkdir(exist_ok=False);(out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
