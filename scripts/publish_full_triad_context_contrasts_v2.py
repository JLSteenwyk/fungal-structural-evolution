#!/usr/bin/env python3
"""Publish all closed context-direction counts and a descriptive reference-policy figure."""
import argparse
import csv
import itertools
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg');matplotlib.rcParams['svg.fonttype']='none'
import matplotlib.pyplot as plt
import numpy as np
from full_triad_context_contrast_sources import MASKS,CORES,POLICIES,PHYSICAL_DIRECTIONS,STATES
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--completion',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    c=json.loads(a.completion.read_text());assert c['status']=='complete_verified_full_triad_context_contrasts' and c['exact_process_journals_checked']==2 and c['scientific_eligibility'] is False
    archive=Path(c['full_hash_archive']);assert sha(archive)==c['full_hash_archive_sha256'];proof=json.loads(archive.read_text())
    assert proof['status']==c['status']+'_archive' and len(proof['services'])==2 and len(proof['source_hashes'])==c['bound_source_hashes']
    rp=Path(c['producer_receipt']);ap=Path(c['independent_readback']);assert sha(rp)==c['producer_receipt_sha256'] and sha(ap)==c['independent_readback_sha256']
    r,v=[json.loads(p.read_text()) for p in [rp,ap]];assert v['producer_receipt_sha256']==sha(rp)
    source=rp.parent/'context_contrast_counts.tsv';assert sha(source)==r['artifacts'][source.name]==proof['source_hashes'][str(source)]
    with source.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==c['summary_rows']==r['summary_rows']==v['summary_rows']==10800
    keys={};screens=set();denominators={}
    for row in rows:
        key=tuple(row[k] for k in ['guide','design','mask_scenario','core_scenario','screen','policy','direction_state']);assert key not in keys
        count=int(row['context_count']);keys[key]=count;screens.add(row['screen']);assert count>=0
        assert count==c['direction_counts'].get('|'.join(key),0)==r['direction_counts'].get('|'.join(key),0)==v['direction_counts'].get('|'.join(key),0)
        baseline=tuple(int(row[k]) for k in ['source_contexts','parent_eligible_contexts','parent_eligible_lexical_gene_present','parent_eligible_lexical_model_present'])
        assert baseline[0]==c['guide_contexts'][row['guide']];pair=key[:2]
        assert pair not in denominators or denominators[pair]==baseline;denominators[pair]=baseline
    assert set(keys)==set(itertools.product(['profile','mafft'],['availability','sequence_first'],MASKS,CORES,screens,POLICIES,STATES)) and len(screens)==6
    for g,d,m,core,s,policy in itertools.product(['profile','mafft'],['availability','sequence_first'],MASKS,CORES,screens,POLICIES):assert sum(keys[g,d,m,core,s,policy,state] for state in STATES)==c['guide_contexts'][g]
    table=Path('docs/tables/full_triad_context_contrast_counts_20261002_v2.tsv');assert not table.exists();table.write_bytes(source.read_bytes())
    pools=POLICIES[2:];states=PHYSICAL_DIRECTIONS+['reference_direction_disagreement'];colors=['#377eb8','#e69f00','#999999','#b85c8a','#333333']
    names=['Positive','Negative','Within numeric tolerance','Sign uncertain','Reference disagreement'];panels=[]
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharey=True);all_numbers=[]
    for ax,(g,d) in zip(axes.ravel(),itertools.product(['profile','mafft'],['availability','sequence_first'])):
        values=np.array([[keys[g,d,'both_masks','both_cores','n50_c70',pool,state] for pool in pools] for state in states]);bottom=np.zeros(3)
        for nums,state,color,label in zip(values,states,colors,names):
            ax.bar(np.arange(3),nums,bottom=bottom,color=color,label=label,width=.6)
            for i,(n,b) in enumerate(zip(nums,bottom)):
                if n and state!='reference_direction_disagreement':ax.text(i,b+n/2,f'{n:,}',ha='center',va='center',fontsize=9,color='#222222' if state=='negative' else 'white')
            bottom+=nums
        for i,(total,disagree) in enumerate(zip(bottom,values[-1])):
            ax.annotate(f'{int(total):,}',xy=(i,total),xytext=(0,5),
                textcoords='offset points',ha='center',fontsize=9)
            if disagree:ax.annotate(f'{int(disagree)} ref. disagreements',
                xy=(i,total),xytext=(0,20),textcoords='offset points',
                ha='center',fontsize=8,color='#333333')
        ax.set_title(g.upper()+' / '+d.replace('_',' '),fontsize=11);ax.set_xticks(np.arange(3),['Fixed lexical','Any eligible ties','All original ties'],fontsize=9);ax.set_ylim(0,6500)
        ax.spines[['top','right']].set_visible(False);panels.append(dict(guide=g,design=d,policies=pools,states=states,counts=values.tolist(),eligible_totals=bottom.astype(int).tolist()))
        all_numbers.extend(int(n) for n in [*bottom,*values.ravel()] if n)
    axes[0,0].set_ylabel('Conditional contexts passing quality');axes[1,0].set_ylabel('Conditional contexts passing quality')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,bbox_to_anchor=(.5,.05),frameon=False,fontsize=9)
    fig.suptitle('Original-context reference-contrast directions\nBoth masks and both cores; ≥50 residues and ≥70% original coverage',fontsize=14)
    fig.text(.5,.015,'Original parent/native/model gates retained. Any-eligible pool agreement need not cover all original ties. Cohorts overlap; no biological significance inferred.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.14,1,.92));outputs=[]
    for extension in ['png','svg','pdf']:
        target=Path('docs/figures/full_triad_context_contrasts_20261002_v2.'+extension);assert not target.exists();fig.savefig(target,dpi=200);outputs.append(target)
    plt.close(fig)
    tree=ET.parse(outputs[1]);ns='{http://www.w3.org/2000/svg}';text=' '.join(''.join(e.itertext()) for e in tree.iter(ns+'text'))
    for n in all_numbers:assert f'{n:,}' in text
    meta=ET.Element(ns+'metadata',id='full-context-contrast-counts');meta.text=json.dumps(panels,separators=(',',':'));tree.getroot().append(meta);tree.write(outputs[1],encoding='utf-8',xml_declaration=True)
    assert json.loads(ET.parse(outputs[1]).find(ns+'metadata[@id="full-context-contrast-counts"]').text)==panels
    result=dict(status='published_full_verified_context_contrast_counts_and_figure',source_completion=str(a.completion),source_completion_sha256=sha(a.completion),full_source_archive=str(archive),full_source_archive_sha256=sha(archive),original_completion_journals=2,
        target_contexts=283409,reference_tie_records=214461,full_summary_rows=10800,plotted_cohorts=4,plotted_policy_bars=12,full_table_partitions_checked=True,svg_text_and_metadata_values_checked=True,
        source_hashes={str(p):sha(p) for p in [a.completion,archive,rp,ap,source,Path(__file__)]},artifacts={str(p):sha(p) for p in [table,*outputs]},scientific_eligibility=False,
        scope='All10800closed disjoint state count rows/283409originalcontexts/214461reference ties, allsource denominators/9mask-core scenarios/6screens/5policies/10states checked. Figure all4guide-design cohorts/3native-both policies underjointmask/core n50_c70; exactcounts/SVGmetadata checked. ActualPNG/PDF visualinspection separate. Overlapping dependent conditional cohorts, not significance, independent events or confidence intervals. Anyeligiblepool agreement not completeoriginaltie agreement. All8aims incomplete; GPU inference paused.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}))


if __name__=='__main__':main()
