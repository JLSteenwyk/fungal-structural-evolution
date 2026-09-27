#!/usr/bin/env python3
"""Verify exact cross-guide identities and all copied fields/directional classifications."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

K=['family','gene_a','gene_b','pfam_accession','screen']
C=['gene_node','taxon_id','chosen_reference_gene','reference_choices','complete','eligible_contrast_min','eligible_contrast_max','duplicate_pair_sequence_distance','sequence_signed_difference','sequence_normalized_contrast','sequence_direction','sequence_covariate_status']

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();p=json.loads(a.plan.read_text());root=Path(p['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
    assert r['status']=='complete_domain_guide_direction_comparison_pending_readback' and r['plan_sha256']==sha(a.plan)
    for path,h in p['pins'].items():assert sha(path)==h,path
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    source=pd.read_csv(Path(p['source'])/'event_domain_robustness.tsv',sep='\t',dtype=str,keep_default_na=False)
    assert not source.duplicated(['guide']+K).any()
    frames=[]
    for guide in ['mafft','profile']:
        frame=source.loc[source.guide.eq(guide),K+C+p['margins']].rename(columns={c:guide+'_'+c for c in C+p['margins']});frames.append(frame)
    joined=frames[0].merge(frames[1],on=K,how='outer',validate='one_to_one',indicator=True)
    for col in joined.columns:
        if col!='_merge':joined[col]=joined[col].fillna('')
    blocks=[]
    for margin in p['margins']:
        e=joined[K+[g+'_'+c for g in ['mafft','profile'] for c in C]].copy();e['margin']=margin
        e['presence']=joined['_merge'].map({'both':'both_guides','left_only':'mafft_only','right_only':'profile_only'}).astype(str)
        both=joined['_merge'].eq('both');struct={};seq={}
        for guide in ['mafft','profile']:
            e[guide+'_structural_direction']=joined[guide+'_'+margin]
            struct[guide]=joined[guide+'_'+margin].map({'a_farther_from_reference':1,'b_farther_from_reference':-1})
            seq[guide]=joined[guide+'_sequence_direction'].map({'a_longer':1,'b_longer':-1}).where(joined[guide+'_sequence_covariate_status'].eq('resolved_pair_distance'))
            present=joined[guide+'_gene_node'].ne('')
            e[guide+'_sequence_structure_relation']=np.select([~present,struct[guide].isna(),seq[guide].isna(),struct[guide].eq(seq[guide])],['missing_guide','structural_direction_unresolved','sequence_direction_unresolved','concordant'],default='discordant')
        e['structural_guide_agreement']=np.select([~both,struct['mafft'].isna()|struct['profile'].isna(),struct['mafft'].eq(struct['profile'])],['missing_guide','unresolved_in_one_or_both','same_stable_direction'],default='opposite_stable_directions')
        e['sequence_guide_agreement']=np.select([~both,seq['mafft'].isna()|seq['profile'].isna(),seq['mafft'].eq(seq['profile'])],['missing_guide','unresolved_in_one_or_both','same_resolved_direction'],default='opposite_resolved_directions')
        x=e['mafft_sequence_structure_relation'];y=e['profile_sequence_structure_relation']
        e['joint_sequence_structure_relation']=np.select([x.eq('concordant')&y.eq('concordant'),x.eq('discordant')&y.eq('discordant'),(x.eq('concordant')&y.eq('discordant'))|(x.eq('discordant')&y.eq('concordant'))],['concordant_in_both_guides','discordant_in_both_guides','guide_sensitive_relationship'],default='unresolved_or_missing')
        blocks.append(e)
    expected=pd.concat(blocks,ignore_index=True).set_index(K+['margin']).sort_index()
    actual=pd.read_csv(root/'guide_comparison.tsv',sep='\t',dtype=str,keep_default_na=False);assert not actual.duplicated(K+['margin']).any();actual=actual.set_index(K+['margin']).sort_index()
    assert expected.index.equals(actual.index) and set(expected.columns)==set(actual.columns)
    for col in expected:assert expected[col].astype(str).equals(actual[col]),col
    flat=actual.reset_index();presence={};counts={}
    for (screen,margin,label),n in flat.groupby(['screen','margin','presence']).size().items():presence['|'.join([screen,margin,label])]=int(n)
    for field in ['structural_guide_agreement','sequence_guide_agreement','mafft_sequence_structure_relation','profile_sequence_structure_relation','joint_sequence_structure_relation']:
        for (screen,margin,label),n in flat.groupby(['screen','margin',field]).size().items():counts['|'.join([screen,margin,field,label])]=int(n)
    assert presence==r['presence_counts'] and counts==r['classification_counts'] and len(actual)==r['comparison_rows'] and len(source)==r['source_rows']
    selected=flat.loc[flat['screen'].eq('n30_c70') & flat['margin'].eq('direction_margin_0_1')]
    stable=selected.loc[selected['structural_guide_agreement'].eq('same_stable_direction')]
    result=dict(status='passed_full_domain_guide_comparison_readback',producer_receipt_sha256=sha(rp),plan_sha256=sha(a.plan),checker_sha256=sha(__file__),rows_checked=len(actual),source_rows_checked=len(source),n30_c70_margin_0_1_presence=selected['presence'].value_counts().to_dict(),n30_c70_margin_0_1_structure_agreement=selected['structural_guide_agreement'].value_counts().to_dict(),stable_structure_sequence_relations=stable['joint_sequence_structure_relation'].value_counts().to_dict(),stable_structural_event_domains=len(stable),distinct_stable_gene_pairs=len(stable[['family','gene_a','gene_b']].drop_duplicates()),scope='Independent outer merge verified all exact gene-pair/Pfam keys, every copied source field and all classification cells/counts across all screens/margins. Directional concordance is descriptive; multiple domains/families/taxa are dependent and no phylogenetic association test was performed.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
