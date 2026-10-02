#!/usr/bin/env python3
"""Publish all closed joint-direction counts and every n50/c70 scenario."""
import argparse
import csv
import itertools
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype']='none'
import matplotlib.pyplot as plt
import numpy as np
from full_triad_joint_direction_sources import MASKS,METHOD_SCENARIOS,CORE_SCENARIOS,DIRECTIONS,QUALIFIED
from screen_duplication_alignment_reuse import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--completion',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    c=json.loads(a.completion.read_text());assert c['status']=='complete_verified_full_triad_joint_directions' and c['scientific_eligibility'] is False and c['exact_process_journals_checked']==2
    archive=Path(c['full_hash_archive']);assert sha(archive)==c['full_hash_archive_sha256'];proof=json.loads(archive.read_text())
    assert proof['status']==c['status']+'_archive' and len(proof['services'])==2 and len(proof['source_hashes'])==c['bound_source_hashes']
    rp,ap=map(Path,[c['producer_receipt'],c['independent_readback']]);assert sha(rp)==c['producer_receipt_sha256'] and sha(ap)==c['independent_readback_sha256']
    r,v=[json.loads(x.read_text()) for x in [rp,ap]];assert v['producer_receipt_sha256']==sha(rp) and r['scientific_eligibility'] is v['scientific_eligibility'] is False
    for k in proof['summary']:assert c[k]==r[k]==v[k]==proof['summary'][k]
    source=rp.parent/'joint_direction_counts.tsv';assert sha(source)==r['artifacts'][source.name]==proof['source_hashes'][str(source)]
    rows=list(csv.DictReader(source.open(),delimiter='\t'));assert len(rows)==c['summary_rows']==1134
    keys={};screens=set()
    for row in rows:
        key=tuple(row[k] for k in ['mask_scenario','method_scenario','core_scenario','screen','qualified_direction']);assert key not in keys
        keys[key]=int(row['physical_triads']);screens.add(row['screen']);assert int(row['total_physical_triads'])==27056 and keys[key]>=0
        for record in [c,r,v]:assert keys[key]==record['qualified_direction_counts'].get('|'.join(key),0)
    assert len(screens)==6 and set(keys)==set(itertools.product(MASKS,METHOD_SCENARIOS,CORE_SCENARIOS,screens,QUALIFIED))
    for m,q,k,s in itertools.product(MASKS,METHOD_SCENARIOS,CORE_SCENARIOS,screens):assert sum(keys[m,q,k,s,d] for d in QUALIFIED)==27056
    table=Path('docs/tables/full_triad_joint_direction_counts_20261002.tsv');assert not table.exists();table.write_bytes(source.read_bytes())
    scenarios=list(itertools.product(MASKS,CORE_SCENARIOS));classes=DIRECTIONS;colors=['#377eb8','#e69f00','#999999','#b85c8a','#6b6b6b','#333333']
    labels=['Positive','Negative','Within numerical tolerance','Sign uncertain','Unavailable','Nonunique fit']
    names=[{'full':'Full','plddt70':'pLDDT ≥70','both_masks':'Both masks'}[m]+'\n'+{'reference_common':'Reference','cycle_consistent':'Cycle','both_cores':'Both cores'}[k] for m,k in scenarios]
    fig,axes=plt.subplots(3,2,figsize=(17,13));all_counts=[];all_passed=[];x=np.arange(9)
    for j,q in enumerate(METHOD_SCENARIOS):
        counts=np.asarray([[keys[m,q,k,'n50_c70',d] for m,k in scenarios] for d in classes]);passed=counts.sum(axis=0);assert (passed>0).all()
        all_counts.append(counts.tolist());all_passed.append(passed.tolist());left,right=axes[j]
        left.bar(x,passed,color='#416b80',width=.75)
        for i,n in enumerate(passed):left.annotate(f'{n:,}',(i,n),xytext=(0,5),textcoords='offset points',ha='center',fontsize=8)
        left.set_ylim(0,max(passed)*1.2);left.set_ylabel('Quality-passing physical triples');left.set_title({'famsa_default':'FAMSA','mafft_auto':'MAFFT','both_methods':'Both sequence methods'}[q]+' + structural correspondence',fontsize=12)
        bottom=np.zeros(9)
        for values,label,color in zip(counts,labels,colors):
            frac=values/passed*100;right.bar(x,frac,bottom=bottom,label=label,color=color,width=.75)
            for i,(n,f,b) in enumerate(zip(values,frac,bottom)):
                if n:right.text(i,b+f/2,f'{n:,}',ha='center',va='center',fontsize=7,color='white' if label!='Negative' else '#222')
            bottom+=frac
        right.set_ylim(0,100);right.set_ylabel('Direction among quality-passing triples (%)');right.set_title('All selected correspondence, mask, core and order alternatives',fontsize=11)
        for ax in [left,right]:ax.set_xticks(x,names,fontsize=8);ax.spines[['top','right']].set_visible(False)
    handles,legend=axes[0,1].get_legend_handles_labels();fig.legend(handles,legend,loc='lower center',bbox_to_anchor=(.5,.035),ncol=3,frameon=False,fontsize=10)
    fig.suptitle('Joint sequence and structural correspondence: reference-distance direction sensitivity',fontsize=17)
    fig.text(.5,.015,'≥50 residues; ≥70% of each original protein. RMSD(A, reference) − RMSD(B, reference). Dependent alternatives; no evolutionary polarity or confidence interval.',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.10,1,.96));outputs=[]
    for ext in ['png','svg','pdf']:
        target=Path('docs/figures/full_triad_joint_directions_20261002.'+ext);assert not target.exists();fig.savefig(target,dpi=160);outputs.append(target)
    plt.close(fig)
    tree=ET.parse(outputs[1]);ns='{http://www.w3.org/2000/svg}';text=' '.join(''.join(e.itertext()) for e in tree.iter(ns+'text'))
    for n in np.asarray(all_counts).ravel().tolist()+np.asarray(all_passed).ravel().tolist():
        if n:assert f'{n:,}' in text
    values=dict(screen='n50_c70',method_scenarios=METHOD_SCENARIOS,mask_core_scenarios=scenarios,directions=classes,counts=all_counts,quality_passing=all_passed)
    item=ET.Element(ns+'metadata',id='joint-direction-counts');item.text=json.dumps(values,separators=(',',':'));tree.getroot().append(item);tree.write(outputs[1],encoding='utf-8',xml_declaration=True)
    assert json.loads(ET.parse(outputs[1]).find(ns+'metadata[@id="joint-direction-counts"]').text)==json.loads(json.dumps(values))
    result=dict(status='published_full_verified_joint_direction_counts_and_figure',source_completion=str(a.completion),source_completion_sha256=sha(a.completion),full_source_archive=str(archive),full_source_archive_sha256=sha(archive),original_completion_journals=2,measured_triads=27056,full_summary_rows=1134,plotted_scenarios=27,figure_screen='n50_c70',plotted_counts_checked_against_full_table=True,svg_text_and_metadata_values_checked=True,
        source_hashes={str(p):sha(p) for p in [a.completion,archive,rp,ap,source,Path(__file__)]},artifacts={str(p):sha(p) for p in [table,*outputs]},scientific_eligibility=False,scope='All1134joint-direction count rows and all27n50c70mask/method/core scenarios published. Exact disjoint denominators/counts checked against original complete producer/reader/two-journal closure; fullsource byte refresh recorded separately. SVG values and metadata verified; actual raster/rendered PDF inspection separate. Source alternatives overlap; direction envelopes notCI, matchingcategories notstable sign or effectagreement, predictedcoordinates sequence-derived, not evolutionarypolarity/calibratedduplicationeffects. All8aims incomplete; GPU predictionpaused.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}),flush=True)


if __name__=='__main__':main()
