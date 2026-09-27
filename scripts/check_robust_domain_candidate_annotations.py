#!/usr/bin/env python3
"""Read back the full annotation join with independent dataframe aggregation."""
import argparse,json,re,hashlib
from pathlib import Path
import pandas as pd
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not a.output.exists()
    p=json.loads(a.plan.read_text());root=Path(p['output']);r=json.loads((root/'receipt.json').read_text());assert r['plan_sha256']==sha(a.plan)
    for f,h in p['pins'].items():assert sha(f)==h
    for f,h in r['artifacts'].items():assert sha(root/f)==h
    source=read(p['candidates']);source=source[(source.screen=='n30_c70')&(source.margin=='direction_margin_0_1')&(source.structural_guide_agreement=='same_stable_direction')]
    keys=['family','gene_a','gene_b','pfam_accession'];out=read(root/'candidate_annotations.tsv');assert len(out)==len(source)==947 and not out.duplicated(keys).any()
    si=source.set_index(keys).sort_index();oi=out.set_index(keys).sort_index();pd.testing.assert_frame_equal(si,oi[si.columns],check_exact=True)
    links=read(p['links']);joined=[]
    for guide in ['mafft','profile']:
        wanted=source[keys+[guide+'_gene_node']].rename(columns={guide+'_gene_node':'gene_node'})
        matched=wanted.merge(links[links.guide==guide],on=keys+['gene_node'],how='left',validate='one_to_many');assert matched.triad_key.notna().all();joined.append(matched[keys+['triad_key']])
    units=pd.concat(joined).drop_duplicates();fit=read(p['fits']);fit=fit[fit.triad_key.isin(units.triad_key)]
    grid=['triad_key','mask','order_ab','order_ar','order_br','mapping_definition'];assert not fit.duplicated(grid).any() and fit.groupby('triad_key').size().eq(32).all()
    assert fit.n30_c70_pass.eq('1').all() and fit.fit_status.eq('computed_unique_at_numeric_tolerance').all()
    metrics=['common_residues','coverage_a','coverage_b','coverage_reference','rmsd_ab','rmsd_ar_minus_br','sequence_identity_ab','mean_plddt_a','mean_plddt_b','mean_plddt_reference','joint_plddt70_fraction']
    numeric=fit[['triad_key']+metrics].copy();numeric[metrics]=numeric[metrics].astype(float);assert np.isfinite(numeric[metrics]).all().all()
    pooled=units.merge(numeric,on='triad_key',validate='many_to_many')
    for metric in metrics:
        for stat in ['min','max']:
            expected=pooled.groupby(keys)[metric].agg(stat).sort_index();assert expected.index.equals(oi.index)
            assert np.allclose(expected.to_numpy(),oi['all_alternatives_'+stat+'_'+metric].astype(float),atol=1e-12,rtol=1e-12)
    sets=units.groupby(keys).triad_key.agg(lambda x:sorted(set(x))).sort_index()
    assert all(json.loads(value)==expected for value,expected in zip(oi.triad_keys_json,sets))
    assert np.array_equal(oi.unique_interval_triads.astype(int),sets.map(len))
    assert np.array_equal(oi.unique_fit_alternatives.astype(int),sets.map(len)*32)
    lo=oi.all_alternatives_min_rmsd_ar_minus_br.astype(float);hi=oi.all_alternatives_max_rmsd_ar_minus_br.astype(float)
    assert ((lo>.1)|(hi<-.1)).all()
    assert np.allclose(oi.minimum_absolute_contrast_angstrom.astype(float),np.minimum(abs(lo),abs(hi)),atol=1e-12,rtol=1e-12)
    meta={}
    for block in Path(p['pfam']).read_text().split('//'):
        values={k:v.strip() for k,v in re.findall(r'^#=GF (AC|ID|DE|TP|CL)\s+(.+)$',block,re.M)}
        if 'AC' in values:assert values['AC'] not in meta;meta[values['AC']]=values
    for _,x in out.iterrows():
        for field,tag in [('pfam_name','ID'),('pfam_description','DE'),('pfam_type','TP'),('pfam_clan','CL')]:assert x[field]==meta[x.pfam_accession].get(tag,'')
    groupkeys=['study_role','candidate_class','pfam_accession'];summary=read(root/'pfam_summary.tsv').set_index(groupkeys);grouped=out.groupby(groupkeys);assert set(summary.index)==set(grouped.groups) and len(summary)==r['pfam_summary_rows']
    for key,xs in grouped:
        row=summary.loc[key];assert int(row.event_domain_combinations)==len(xs)
        assert int(row.distinct_gene_pairs)==len(xs[['family','gene_a','gene_b']].drop_duplicates()) and int(row.families)==xs.family.nunique() and int(row.taxa)==xs.taxon.nunique()
        assert row.pfam_name==meta[key[2]]['ID'] and row.pfam_description==meta[key[2]]['DE']
    assert fit.triad_key.nunique()==r['unique_interval_triads'] and len(fit)==r['unique_fit_alternatives']
    discordant=out[(out.candidate_class=='stable_discordant')&(out.study_role=='ingroup')]
    overview=[]
    for pfam,xs in discordant.groupby('pfam_accession'):
        overview.append(dict(pfam_accession=pfam,name=meta[pfam]['ID'],description=meta[pfam]['DE'],event_domain_combinations=len(xs),families=xs.family.nunique(),taxa=xs.taxon.nunique()))
    result=dict(status='passed_full_robust_domain_candidate_annotation_readback',candidates=len(out),unique_interval_triads=fit.triad_key.nunique(),unique_fit_alternatives=len(fit),pfam_summary_rows=len(summary),plan_sha256=sha(a.plan),producer_receipt_sha256=sha(root/'receipt.json'),checker_sha256=sha(__file__),fungal_discordant_domain_counts=sorted(overview,key=lambda x:(-x['event_domain_combinations'],x['pfam_accession'])),scope='All copied fields, Pfam labels, exact triad sets, fit extrema and domain summary denominators independently joined/aggregated. Descriptive candidate annotations, not independent events or validated functions.')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='fungal_discordant_domain_counts'},indent=2));print('Most represented fungal discordant domains:',result['fungal_discordant_domain_counts'][:5])


if __name__=='__main__':main()
