#!/usr/bin/env python3
"""Publish every closed context joint-direction count and all method/reference cohorts."""
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
from full_triad_context_joint_direction_sources import MASKS,METHOD_SCENARIOS,CORE_SCENARIOS,POLICIES,PHYSICAL_DIRECTIONS,STATES
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--completion',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    c=json.loads(a.completion.read_text());assert c['status']=='complete_verified_full_triad_context_joint_directions' and c['exact_process_journals_checked']==2 and c['scientific_eligibility'] is False
    archive=Path(c['full_hash_archive']);assert sha(archive)==c['full_hash_archive_sha256'];proof=json.loads(archive.read_text())
    assert proof['status']==c['status']+'_archive' and len(proof['services'])==2 and len(proof['source_hashes'])==c['bound_source_hashes']
    rp=Path(c['producer_receipt']);ap=Path(c['independent_readback']);assert sha(rp)==c['producer_receipt_sha256'] and sha(ap)==c['independent_readback_sha256']
    r,v=[json.loads(p.read_text()) for p in [rp,ap]];assert v['producer_receipt_sha256']==sha(rp)
    source=rp.parent/'context_joint_direction_counts.tsv';assert sha(source)==r['artifacts'][source.name]==proof['source_hashes'][str(source)]
    with source.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==c['summary_rows']==r['summary_rows']==v['summary_rows']==32400
    keys={};screens=set();denominators={}
    for row in rows:
        key=tuple(row[k] for k in ['guide','design','mask_scenario','method_scenario','core_scenario','screen','policy','direction_state']);assert key not in keys
        count=int(row['context_count']);keys[key]=count;screens.add(row['screen']);assert count>=0
        assert count==c['direction_counts'].get('|'.join(key),0)==r['direction_counts'].get('|'.join(key),0)==v['direction_counts'].get('|'.join(key),0)
        baseline=tuple(int(row[k]) for k in ['source_contexts','parent_eligible_contexts','parent_eligible_lexical_gene_present','parent_eligible_lexical_model_present'])
        assert baseline[0]==c['guide_contexts'][row['guide']];pair=key[:2]
        assert pair not in denominators or denominators[pair]==baseline;denominators[pair]=baseline
    assert set(keys)==set(itertools.product(['profile','mafft'],['availability','sequence_first'],MASKS,METHOD_SCENARIOS,CORE_SCENARIOS,screens,POLICIES,STATES)) and len(screens)==6
    for g,d,m,method,core,s,policy in itertools.product(['profile','mafft'],['availability','sequence_first'],MASKS,METHOD_SCENARIOS,CORE_SCENARIOS,screens,POLICIES):assert sum(keys[g,d,m,method,core,s,policy,state] for state in STATES)==c['guide_contexts'][g]
    table=Path('docs/tables/full_triad_context_joint_direction_counts_20261002.tsv');assert not table.exists();table.write_bytes(source.read_bytes())
    pools=POLICIES[2:];states=PHYSICAL_DIRECTIONS+['reference_direction_disagreement'];colors=['#377eb8','#e69f00','#999999','#b85c8a','#333333']
    names=['Positive','Negative','Within numeric tolerance','Sign uncertain','Reference disagreement'];panels=[]
    fig,axes=plt.subplots(4,3,figsize=(16,15),sharey=True);all_numbers=[];maximum=0
    cohorts=list(itertools.product(['profile','mafft'],['availability','sequence_first']))
    for row,(g,d) in enumerate(cohorts):
        for column,method in enumerate(METHOD_SCENARIOS):
            ax=axes[row,column]
            values=np.array([[keys[g,d,'both_masks',method,'both_cores','n50_c70',pool,state] for pool in pools] for state in states]);bottom=np.zeros(3)
            for nums,state,color,label in zip(values,states,colors,names):
                ax.bar(np.arange(3),nums,bottom=bottom,color=color,label=label,width=.6)
                for i,(n,b) in enumerate(zip(nums,bottom)):
                    if n and state!='reference_direction_disagreement':ax.text(i,b+n/2,f'{n:,}',ha='center',va='center',fontsize=9,color='#222222' if state=='negative' else 'white')
                bottom+=nums
            for i,(total,disagree) in enumerate(zip(bottom,values[-1])):
                ax.annotate(f'{int(total):,}',xy=(i,total),xytext=(0,5),textcoords='offset points',ha='center',fontsize=9)
                if disagree:ax.annotate(f'{int(disagree)} ref. disagreements',xy=(i,total),xytext=(0,21),textcoords='offset points',ha='center',fontsize=8,color='#333333')
            maximum=max(maximum,float(max(bottom)))
            ax.set_title(g.upper()+' / '+d.replace('_',' ')+' / '+{'famsa_default':'FAMSA','mafft_auto':'MAFFT','both_methods':'Both methods'}[method],fontsize=11)
            ax.set_xticks(np.arange(3),['Fixed lexical','Any eligible ties','All original ties'],fontsize=9)
            ax.spines[['top','right']].set_visible(False)
            panels.append(dict(guide=g,design=d,method_scenario=method,policies=pools,states=states,counts=values.tolist(),eligible_totals=bottom.astype(int).tolist()))
            all_numbers.extend(int(n) for n in [*bottom,*values.ravel()] if n)
    for ax in axes.ravel():ax.set_ylim(0,max(100,maximum*1.2))
    for ax in axes[:,0]:ax.set_ylabel('Conditional contexts passing quality')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,bbox_to_anchor=(.5,.043),frameon=False,fontsize=10)
    fig.suptitle('Joint sequence and structural correspondences in original gene contexts\nBoth masks and both cores; ≥50 residues and ≥70% original coverage',fontsize=16)
    fig.text(.5,.015,'All original source gates retained. Eligible-pool agreement may omit original ties. Cohorts overlap; directions are not evolutionary polarity or significance.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.11,1,.94));outputs=[]
    for extension in ['png','svg','pdf']:
        target=Path('docs/figures/full_triad_context_joint_directions_20261002.'+extension);assert not target.exists();fig.savefig(target,dpi=160);outputs.append(target)
    plt.close(fig)
    tree=ET.parse(outputs[1]);ns='{http://www.w3.org/2000/svg}';text=' '.join(''.join(e.itertext()) for e in tree.iter(ns+'text'))
    for n in all_numbers:assert f'{n:,}' in text
    meta=ET.Element(ns+'metadata',id='full-context-joint-direction-counts');meta.text=json.dumps(panels,separators=(',',':'));tree.getroot().append(meta);tree.write(outputs[1],encoding='utf-8',xml_declaration=True)
    assert json.loads(ET.parse(outputs[1]).find(ns+'metadata[@id="full-context-joint-direction-counts"]').text)==panels
    result=dict(status='published_full_verified_context_joint_direction_counts_and_figure',source_completion=str(a.completion),source_completion_sha256=sha(a.completion),full_source_archive=str(archive),full_source_archive_sha256=sha(archive),original_completion_journals=2,
        target_contexts=283409,reference_tie_records=214461,full_summary_rows=32400,plotted_cohorts=4,plotted_method_scenarios=3,plotted_policy_bars=36,full_table_partitions_checked=True,svg_text_and_metadata_values_checked=True,
        source_hashes={str(p):sha(p) for p in [a.completion,archive,rp,ap,source,Path(__file__)]},artifacts={str(p):sha(p) for p in [table,*outputs]},scientific_eligibility=False,
        scope='All32400closed disjoint state count rows/283409originalcontexts/214461reference ties, allsource denominators/27mask-method-core scenarios/6screens/5policies/10states checked. Figure all4guide-design cohorts/3sequence-method alternatives/3native-both policies underjointmask/core n50_c70; exactcounts/SVGmetadata checked. ActualPNG/PDF visualinspection separate. Overlapping dependent conditional cohorts, not significance, independent events or confidence intervals. Anyeligiblepool agreement not completeoriginaltie agreement. All8aims incomplete; GPU inference paused.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}))


if __name__=='__main__':main()
