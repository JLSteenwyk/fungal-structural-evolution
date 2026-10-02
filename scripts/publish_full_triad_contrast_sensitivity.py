#!/usr/bin/env python3
"""Publish complete closed contrast counts and a descriptive sensitivity figure."""
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
from full_triad_contrast_sensitivity_sources import MASKS,CORES,QUALIFIED
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--completion',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    completed=json.loads(a.completion.read_text())
    assert completed['status']=='complete_verified_full_triad_contrast_sensitivity' and completed['exact_process_journals_checked']==2
    assert completed['scientific_eligibility'] is False
    archive=Path(completed['full_hash_archive']);assert sha(archive)==completed['full_hash_archive_sha256'];proof=json.loads(archive.read_text())
    assert proof['status']==completed['status']+'_archive' and len(proof['services'])==2 and len(proof['source_hashes'])==completed['bound_source_hashes']
    rp=Path(completed['producer_receipt']);reader=Path(completed['independent_readback'])
    assert sha(rp)==completed['producer_receipt_sha256'] and sha(reader)==completed['independent_readback_sha256']
    r=json.loads(rp.read_text());checked=json.loads(reader.read_text());assert checked['producer_receipt_sha256']==sha(rp)
    assert r['scientific_eligibility'] is checked['scientific_eligibility'] is False
    source=rp.parent/'contrast_sensitivity_counts.tsv';assert sha(source)==r['artifacts'][source.name]==proof['source_hashes'][str(source)]
    with source.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==completed['summary_rows']==r['summary_rows']==checked['summary_rows']==378
    keys={};screens=set()
    for row in rows:
        key=tuple(row[k] for k in ['mask_scenario','core_scenario','screen','qualified_direction'])
        assert key not in keys;keys[key]=int(row['physical_triads']);screens.add(row['screen'])
        assert int(row['total_physical_triads'])==27056 and keys[key]>=0
        assert keys[key]==completed['qualified_direction_counts'].get('|'.join(key),0)==r['qualified_direction_counts'].get('|'.join(key),0)==checked['qualified_direction_counts'].get('|'.join(key),0)
    assert len(screens)==6 and set(keys)==set(itertools.product(MASKS,CORES,screens,QUALIFIED))
    for m,c,s in itertools.product(MASKS,CORES,screens):assert sum(keys[m,c,s,d] for d in QUALIFIED)==27056
    table=Path('docs/tables/full_triad_contrast_sensitivity_counts_20261002.tsv');assert not table.exists();table.write_bytes(source.read_bytes())
    scenarios=list(itertools.product(MASKS,CORES));classes=['positive','negative','sign_uncertain','within_numerical_tolerance','unavailable','nonunique_fit']
    counts=np.array([[keys[m,c,'n50_c70',d] for m,c in scenarios] for d in classes]);passed=counts.sum(axis=0);assert all(passed>0)
    colors=['#377eb8','#e69f00','#b85c8a','#999999','#6b6b6b','#333333']
    labels=['Positive','Negative','Sign uncertain','Within numeric tolerance','Unavailable','Nonunique fit']
    names=[f'{m.replace("plddt70","pLDDT ≥70").replace("both_masks","Both masks").replace("full","Full")}\n{c.replace("reference_common","Reference common").replace("cycle_consistent","Cycle consistent").replace("both_cores","Both cores")}' for m,c in scenarios]
    fig,axes=plt.subplots(2,1,figsize=(12,8),sharex=True,gridspec_kw={'height_ratios':[1,1]});x=np.arange(9)
    axes[0].bar(x,passed,color='#416b80',width=.7)
    for i,n in enumerate(passed):axes[0].text(i,n+110,f'{n:,}',ha='center',fontsize=9)
    axes[0].set_ylim(0,passed.max()*1.2);axes[0].set_ylabel('Physical triples passing quality');axes[0].set_title('All structural orders: ≥50 residues and ≥70% of each original protein',fontsize=12)
    bottom=np.zeros(9)
    for values,label,color in zip(counts,labels,colors):
        frac=values/passed*100;axes[1].bar(x,frac,bottom=bottom,label=label,color=color,width=.7)
        for i,(n,pct,b) in enumerate(zip(values,frac,bottom)):
            if n:axes[1].text(i,b+pct/2,f'{n:,}',ha='center',va='center',fontsize=8,color='white' if label!='Negative' else '#222222')
        bottom+=frac
    axes[1].set_ylim(0,100);axes[1].set_ylabel('Direction among quality-passing triples (%)');axes[1].set_xticks(x,names,fontsize=8)
    axes[1].legend(loc='upper center',bbox_to_anchor=(.5,-.20),ncol=3,frameon=False,fontsize=9)
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Reference-distance contrast sensitivity across masks and cores',fontsize=15)
    fig.text(.5,.015,'RMSD(A, reference) − RMSD(B, reference); numerical boundary 10⁻⁹ Å. Dependent alternatives, not evolutionary effects.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.075,1,.96));outputs=[]
    for extension in ['png','svg','pdf']:
        target=Path('docs/figures/full_triad_contrast_sensitivity_20261002.'+extension);assert not target.exists();fig.savefig(target,dpi=200);outputs.append(target)
    plt.close(fig)
    svg=outputs[1];tree=ET.parse(svg);ns='{http://www.w3.org/2000/svg}'
    text=' '.join(''.join(e.itertext()) for e in tree.iter(ns+'text'))
    for n in [*passed,*counts.ravel()]:
        if n:assert f'{n:,}' in text
    metadata=ET.Element(ns+'metadata',id='full-contrast-counts');metadata.text=json.dumps(dict(screen='n50_c70',scenarios=scenarios,classes=classes,counts=counts.tolist(),quality_passing=passed.tolist()),separators=(',',':'))
    tree.getroot().append(metadata);tree.write(svg,encoding='utf-8',xml_declaration=True)
    parsed=ET.parse(svg).find(ns+'metadata[@id="full-contrast-counts"]');assert json.loads(parsed.text)['counts']==counts.tolist()
    bound=[a.completion,archive,rp,reader,source,Path(__file__)]
    result=dict(status='published_full_verified_structural_contrast_sensitivity_counts_and_figure',source_completion=str(a.completion),source_completion_sha256=sha(a.completion),
        full_source_archive=str(archive),full_source_archive_sha256=sha(archive),original_completion_journals=2,measured_triads=27056,full_summary_rows=378,
        plotted_scenarios=9,figure_screen='n50_c70',plotted_counts_checked_against_full_table=True,svg_text_and_metadata_values_checked=True,
        source_hashes={str(p):sha(p) for p in bound},artifacts={str(p):sha(p) for p in [table,*outputs]},scientific_eligibility=False,
        scope='Full378row/9scenario/6screen/7category disjoint descriptive physical contrast counts. Plotall9n50_c70scenarios, exactqualitypassed counts and allnonzero qualifieddirectioncounts; SVG text/metadata checked. Publication requires complete original two-journal source closure. Raster/PDF visual inspection recorded separately. Sensitivity alternatives are dependent, envelopes are not confidence intervals, reference-distance direction is not evolutionary polarity/asymmetry significance. All8aims incomplete; GPU prediction paused.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}))


if __name__=='__main__':main()
