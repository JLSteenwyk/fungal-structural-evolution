#!/usr/bin/env python3
"""Independently check complete triad joins and categories with dataframe merges."""
import argparse,json,hashlib
from pathlib import Path
import pandas as pd


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    recpath=a.source/'receipt.json';r=json.loads(recpath.read_text())
    if r['status']!='complete_duplication_triad_architecture_controls':raise ValueError('Wrong producer')
    for path,h in r['source_sha256'].items():
        if sha(path)!=h:raise ValueError('Changed source')
    for name,h in r['artifacts'].items():
        if sha(a.source/name)!=h:raise ValueError('Changed output')
    def source(name):
        candidates=[Path(x) for x in r['source_sha256'] if Path(x).name==name]
        if len(candidates)!=1:raise ValueError('Ambiguous source')
        return pd.read_csv(candidates[0],sep='\t',dtype=str,keep_default_na=False)
    ref=source('event_reference_comparisons.tsv');base=source('event_model_pair_links.tsv');controls=source('pair_architecture_controls.tsv')
    observed=pd.read_csv(a.source/'triad_architecture.tsv',sep='\t',dtype=str,keep_default_na=False)
    event=['guide','family','gene_node','gene_a','gene_b'];triad=event+['reference_gene'];keys=triad+['policy']
    policies=pd.DataFrame({'policy':['alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore']})
    expected=ref[ref.focal_side.eq('a')][triad].merge(policies,how='cross')
    if observed.duplicated(keys).any() or expected.duplicated(keys).any():raise ValueError('Duplicate triad/policy')
    scope=expected.merge(observed[keys],on=keys,how='outer',indicator=True,validate='one_to_one')
    if not scope._merge.eq('both').all():raise ValueError('Incomplete triad/policy grid')
    for side in ['a','b']:
        src=ref[ref.focal_side.eq(side)]
        check=observed.merge(src,on=triad,suffixes=('','_source'),validate='many_to_one')
        mappings={f'model_{side}':'focal_model',f'version_{side}':'focal_version',f'{side}_reference_pair_key':'pair_key',
                  'reference_model':'reference_model_source','reference_version':'reference_version_source','lexical_representative':'lexical_representative_source'}
        for k,v in mappings.items():
            if not check[k].eq(check[v]).all():raise ValueError('Reference join mismatch: '+k)
    check=observed.merge(base,on=event,suffixes=('','_source'),validate='many_to_one')
    for k,v in [('duplicate_pair_key','pair_key'),('model_a','model_a_source'),('model_b','model_b_source'),('version_a','version_a_source'),('version_b','version_b_source')]:
        if not check[k].eq(check[v]).all():raise ValueError('Duplicate join mismatch: '+k)
    identities=[]
    for label,left,right in [('duplicate','a','b'),('a_reference','a','reference'),('b_reference','b','reference')]:
        def col(prefix):return (observed['reference_model'],observed['reference_version']) if prefix=='reference' else (observed['model_'+prefix],observed['version_'+prefix])
        lm,lv=col(left);rm,rv=col(right);same=lm.eq(rm)&lv.eq(rv);identities.append(same)
        if not observed.loc[same,label+'_annotation_class'].eq('identical_model').all() or not observed.loc[same,label+'_both_conservative'].eq('').all():raise ValueError('Identity disposition mismatch')
        check=observed.loc[~same].merge(controls,left_on=[label+'_pair_key','policy'],right_on=['pair_key','policy'],how='left',validate='many_to_one',indicator=True)
        if not check._merge.eq('both').all():raise ValueError('Missing pair control')
        for out,src in [(label+'_annotation_class','annotation_class'),(label+'_both_conservative','both_conservative')]:
            if not check[out].eq(check[src]).all():raise ValueError('Pair control join mismatch')
    cols=[x+'_annotation_class' for x in ['duplicate','a_reference','b_reference']]
    flags=[x+'_both_conservative' for x in ['duplicate','a_reference','b_reference']]
    classification=pd.Series('conservative_same_ordered',index=observed.index)
    classification.loc[~observed[flags].eq('1').all(axis=1)]='same_ordered_not_conservative'
    classification.loc[observed[cols].isin(['different_annotation_content','same_content_different_order']).any(axis=1)]='annotation_difference'
    classification.loc[observed[cols].isin(['neither_annotated','one_unannotated']).any(axis=1)]='incomplete_annotation'
    classification.loc[identities[0]|identities[1]|identities[2]]='identical_model_in_triad'
    if not classification.eq(observed.triad_class).all():raise ValueError('Triad category mismatch')
    counts={':'.join(k):int(v) for k,v in observed.groupby(['guide','policy','triad_class']).size().items()}
    if counts!=r['counts']:raise ValueError('Count mismatch')
    stable=observed.groupby(triad).triad_class.nunique()
    stability={'same_class_all_policies':int(stable.eq(1).sum()),'policy_sensitive_class':int(stable.gt(1).sum())}
    if stability!=r['policy_stability'] or len(stable)!=r['triads'] or len(observed)!=r['policy_rows']:raise ValueError('Scope/stability mismatch')
    result=dict(status='passed_full_duplication_triad_architecture_readback',producer_receipt_sha256=sha(recpath),triads=len(stable),policy_rows=len(observed),policy_stability=stability,script_sha256=sha(__file__),scope='Independent complete-grid dataframe joins to source event/reference ledgers and audited pair controls, identity dispositions, classification and totals. Does not rerun Pfam annotation or infer biological events.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
