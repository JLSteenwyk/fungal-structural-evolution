#!/usr/bin/env python3
"""Annotate all guide comparisons and expose lineage/taxon/family concentration without pooling events."""
import argparse,hashlib,json
from pathlib import Path
import pandas as pd


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    a=ap.parse_args();p=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed plan')
        for path,h in p['pins'].items():
            if sha(path)!=h:raise ValueError('Changed source: '+path)
    verify();source=Path(p['comparison']);sr=json.loads((source/'receipt.json').read_text());audit=json.loads(Path(p['readback']).read_text())
    if audit['status']!='passed_full_domain_guide_comparison_readback' or audit['producer_receipt_sha256']!=sha(source/'receipt.json'):raise ValueError('Unbound comparison audit')
    table=source/'guide_comparison.tsv'
    if sha(table)!=sr['artifacts'][table.name]:raise ValueError('Changed comparison')
    membership=Path(p['membership']);mr=json.loads((membership/'receipt.json').read_text());mt=membership/'taxon_coverage_with_membership.tsv'
    if mr['status']!='complete_duplication_coverage_reconciliation_membership' or sha(mt)!=mr['artifacts'][mt.name]:raise ValueError('Changed membership')
    members=pd.read_csv(mt,sep='\t',dtype=str,keep_default_na=False)
    tax=members.loc[members.in_reconciliation.eq('1'),['taxon','species_name','study_role','lineage']].drop_duplicates()
    if tax.taxon.duplicated().any() or len(tax)!=526:raise ValueError('Inconsistent reconciled taxonomy')
    tax['lineage_group']=tax.lineage.str.split(';').str[0]
    data=pd.read_csv(table,sep='\t',dtype=str,keep_default_na=False)
    data['taxon']=data.mafft_taxon_id.where(data.mafft_taxon_id.ne(''),data.profile_taxon_id)
    annotated=data.merge(tax,on='taxon',how='left',validate='many_to_one')
    if annotated.species_name.isna().any() or len(annotated)!=sr['comparison_rows']:raise ValueError('Unmapped candidate taxon')
    annotated['candidate_class']='structural_direction_unresolved'
    annotated.loc[annotated.presence.ne('both_guides'),'candidate_class']='guide_specific_domain'
    stable=annotated.structural_guide_agreement.eq('same_stable_direction')
    annotated.loc[stable,'candidate_class']='stable_structure_sequence_unresolved'
    annotated.loc[stable & annotated.joint_sequence_structure_relation.eq('concordant_in_both_guides'),'candidate_class']='stable_concordant'
    annotated.loc[stable & annotated.joint_sequence_structure_relation.eq('discordant_in_both_guides'),'candidate_class']='stable_discordant'
    annotated.loc[stable & annotated.joint_sequence_structure_relation.eq('guide_sensitive_relationship'),'candidate_class']='stable_structure_sequence_guide_sensitive'
    annotated['gene_pair_key']=annotated[['family','gene_a','gene_b']].apply(lambda row:json.dumps(row.tolist(),separators=(',',':')),axis=1)
    summaries=[];family=[]
    classes=['structural_direction_unresolved','guide_specific_domain','stable_structure_sequence_unresolved','stable_concordant','stable_discordant','stable_structure_sequence_guide_sensitive']
    base=tax.groupby(['study_role','lineage_group']).agg(reconciled_taxa=('taxon','nunique')).reset_index()
    for (screen,margin),cohort in annotated.groupby(['screen','margin'],sort=True):
        for category in classes:
            selected=cohort.loc[cohort.candidate_class.eq(category)]
            group=selected.groupby(['study_role','lineage_group']).agg(event_domain_combinations=('pfam_accession','size'),distinct_gene_pairs=('gene_pair_key','nunique'),families=('family','nunique'),pfams=('pfam_accession','nunique'),taxa_represented=('taxon','nunique')).reset_index()
            group=base.merge(group,on=['study_role','lineage_group'],how='left',validate='one_to_one').fillna(0)
            for c in ['event_domain_combinations','distinct_gene_pairs','families','pfams','taxa_represented']:group[c]=group[c].astype(int)
            group['screen']=screen;group['margin']=margin;group['candidate_class']=category;summaries.append(group)
            fam=selected.groupby('family').agg(event_domain_combinations=('pfam_accession','size'),distinct_gene_pairs=('gene_pair_key','nunique'),taxa_represented=('taxon','nunique'),pfams=('pfam_accession','nunique')).reset_index()
            fam['screen']=screen;fam['margin']=margin;fam['candidate_class']=category;family.append(fam)
    out=Path(p['output']);out.mkdir(parents=True,exist_ok=False)
    annotated.to_csv(out/'annotated_comparisons.tsv',sep='\t',index=False)
    lineages=pd.concat(summaries,ignore_index=True);lineages.to_csv(out/'lineage_summary.tsv',sep='\t',index=False)
    families=pd.concat(family,ignore_index=True);families.to_csv(out/'family_summary.tsv',sep='\t',index=False)
    selected=annotated.loc[annotated.screen.eq('n30_c70') & annotated.margin.eq('direction_margin_0_1') & stable]
    # Re-evaluate using the annotated table index; no alternate sample selection.
    overview=[]
    for category,part in selected.groupby('candidate_class'):
        overview.append(dict(candidate_class=category,event_domain_combinations=len(part),distinct_gene_pairs=part.gene_pair_key.nunique(),families=part.family.nunique(),pfams=part.pfam_accession.nunique(),taxa=part.taxon.nunique(),manifest_lineage_groups=part.lineage_group.nunique(),study_roles=part.study_role.value_counts().to_dict()))
    verify();result=dict(status='complete_domain_candidate_sampling_pending_readback',plan_sha256=ph,annotated_rows=len(annotated),lineage_summary_rows=len(lineages),family_summary_rows=len(families),reconciled_taxa=len(tax),n30_c70_margin_0_1_overview=overview,artifacts={f.name:sha(f) for f in out.iterdir() if f.is_file()},scope='All screens/margins and guide-comparison rows retained. Manifest broad-lineage annotations and family/taxon/domain multiplicity counts; reconciled taxa form an availability denominator, not a sampling-bias correction. Zero rows explicit. No enrichment, phylogenetic association or fungal-wide generalization.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
