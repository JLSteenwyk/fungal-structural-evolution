#!/usr/bin/env python3
"""Recompute taxon/family coverage summaries from every exported event disposition."""
import argparse,json,math
from pathlib import Path
import pandas as pd
from run_ortholog_pair_guide_comparison import sha

COLUMNS=['terminal_singleton_events','neither_model','one_model','both_models','both_same_model']


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plan=json.loads(a.source_plan.read_text());ph=sha(a.source_plan);folder=Path(plan['output']);receipt=json.loads((folder/'receipt.json').read_text());rh=sha(folder/'receipt.json')
    def verify():
        if sha(a.source_plan)!=ph or sha(folder/'receipt.json')!=rh:raise ValueError('Changed source plan/receipt')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed input: '+path)
        for name,h in receipt['artifacts'].items():
            if sha(folder/name)!=h:raise ValueError('Changed output artifact')
    verify()
    if receipt['status']!='complete_terminal_duplication_sampling_coverage' or receipt['plan_sha256']!=ph:raise ValueError('Unbound coverage receipt')
    read=lambda path:pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False)
    events=read(folder/'terminal_event_coverage.tsv.gz');taxa=read(plan['sampling']).set_index('taxon_id',verify_integrity=True)
    if events.duplicated(['guide','family','taxon_id','gene_node']).any():raise ValueError('Repeated terminal event identity')
    if not set(events.guide)=={'profile','mafft'} or not set(events.taxon_id)<=set(taxa.index):raise ValueError('Unknown event guide/taxon')
    # Reclassify coverage directly from exported model identities, not stored labels.
    n=events.model_a.ne('').astype(int)+events.model_b.ne('').astype(int)
    labels=n.map({0:'neither_model',1:'one_model',2:'both_models'})
    same=events.model_a.ne('') & events.model_a.eq(events.model_b)
    if not labels.equals(events.coverage_class) or not same.astype(int).astype(str).equals(events.same_model):raise ValueError('Event disposition/model identity differs')
    events['terminal_singleton_events']=1
    for label in ['neither_model','one_model','both_models']:events[label]=labels.eq(label).astype(int)
    events['both_same_model']=same.astype(int)
    checked={}
    for level,idcol in [('taxon','taxon_id'),('family','family')]:
        expected=events.groupby(['guide',idcol],sort=True)[COLUMNS].sum()
        if level=='taxon':
            expected=expected.reindex(pd.MultiIndex.from_product([['profile','mafft'],taxa.index],names=['guide',idcol]),fill_value=0)
        output=read(folder/(level+'_coverage.tsv')).rename(columns={level:idcol})
        if output.duplicated(['guide',idcol]).any():raise ValueError('Repeated summary group')
        output=output.set_index(['guide',idcol])
        if set(output.index)!=set(expected.index):raise ValueError('Missing/extra summary group')
        output=output.loc[expected.index]
        if not (output[COLUMNS].astype(int).to_numpy()==expected.to_numpy()).all():raise ValueError('Summary count differs')
        for index,row in output.iterrows():
            total=int(row.terminal_singleton_events);numerator=int(row.both_models)
            if not total:
                if row.both_model_fraction!='':raise ValueError('Zero denominator fraction not blank')
            elif not math.isclose(float(row.both_model_fraction),numerator/total,rel_tol=0,abs_tol=1e-15):raise ValueError('Coverage fraction differs')
            if level=='taxon' and any(row[k]!=taxa.loc[index[1],k] for k in ['species_name','study_role','lineage']):raise ValueError('Taxon metadata differs')
        checked[level]=len(output)
    summaries=[]
    for guide,part in events.groupby('guide'):
        old=next(g for g in receipt['guides'] if g['guide']==guide);totals={k:int(part[k].sum()) for k in COLUMNS}
        if totals!=old['counts']:raise ValueError('Receipt totals differ')
        summary=dict(guide=guide,counts=totals,taxa_with_terminal_singletons=part.taxon_id.nunique(),taxa_with_both_models=part.loc[part.both_models.eq(1),'taxon_id'].nunique(),families_with_terminal_singletons=part.family.nunique(),families_with_both_models=part.loc[part.both_models.eq(1),'family'].nunique())
        if any(old[k]!=v for k,v in summary.items()):raise ValueError('Receipt group universe differs')
        summaries.append(summary)
    verify();result=dict(status='passed_complete_duplication_sampling_aggregate_readback',producer_receipt_sha256=rh,event_rows=len(events),summary_rows=checked,guides=summaries,scope='Every exported model-presence/same-model disposition reclassified and all taxon/family counts, fractions, metadata and summary universes independently recomputed using grouped event data. Native event selection and original bridge join were checked by the producer, not independently repeated here. No missingness correction, phylogenetic model or biological duplication validation.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(event_rows=len(events),summary_rows=checked)),flush=True)


if __name__=='__main__':main()
