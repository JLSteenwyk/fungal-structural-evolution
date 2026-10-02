#!/usr/bin/env python3
"""Publish every fixed-matching balance group with all scenario extrema retained."""
import argparse
import csv
from collections import defaultdict
import itertools
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype']='none'
import matplotlib.pyplot as plt
import numpy as np
from screen_duplication_alignment_reuse import sha

METRICS=['standardized_mean_difference','selection_shift_in_baseline_sd','metadata_matched_shift_in_baseline_sd','target_eligible_all_shift_in_baseline_sd']


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--completion',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    c=json.loads(a.completion.read_text());assert c['status']=='complete_verified_full_fixed_matching_post_screen_balance_and_reuse' and c['exact_process_journals_checked']==2 and c['scientific_eligibility'] is False
    archive=Path(c['full_hash_archive']);assert sha(archive)==c['full_hash_archive_sha256'];proof=json.loads(archive.read_text());assert proof['status']=='complete_verified_full_screened_balance_archive' and len(proof['services'])==2 and len(proof['source_hashes'])==c['bound_source_hashes']
    rp,ap=map(Path,[c['producer_receipt'],c['independent_readback']]);assert sha(rp)==c['producer_receipt_sha256'] and sha(ap)==c['independent_readback_sha256'];r,v=[json.loads(p.read_text()) for p in [rp,ap]]
    assert v['producer_receipt_sha256']==sha(rp) and r['scientific_eligibility'] is v['scientific_eligibility'] is False
    for k,val in proof['summary'].items():assert val==c[k]==r[k]==v[k]
    paths=[rp.parent/'balance.tsv',rp.parent/'coverage.tsv']
    for q in paths:assert sha(q)==r['artifacts'][q.name]==proof['source_hashes'][str(q)]
    rows=list(csv.DictReader(paths[0].open(),delimiter='\t'));coverage=list(csv.DictReader(paths[1].open(),delimiter='\t'))
    assert len(rows)==c['balance_rows']==62208 and len(coverage)==c['coverage_rows']==7776
    guides=c['guides'];policies=c['policies'];masks=c['masks'];screens=[s['id'] for s in c['screens']];features=c['features'];scenario_ids=sorted({row['scenario_id'] for row in rows});assert len(scenario_ids)==54
    coverage_keys={};source_keys=set();groups=defaultdict(list)
    for row in coverage:
        key=tuple(row[k] for k in ['guide','policy','scenario_id','mask','screen']);assert key not in coverage_keys;coverage_keys[key]=row
        assert int(row['metadata_matched_records'])+int(row['metadata_unmatched_records'])==int(row['all_target_records'])
        assert int(row['retained_matches'])+int(row['screen_excluded_matches'])==int(row['metadata_matched_records'])
    assert set(coverage_keys)==set(itertools.product(guides,policies,scenario_ids,masks,screens))
    for row in rows:
        key=tuple(row[k] for k in ['guide','policy','scenario_id','mask','screen','feature']);assert key not in source_keys;source_keys.add(key)
        assert int(row['pairs'])+int(row['retained_feature_excluded_pairs'])==int(coverage_keys[key[:5]]['retained_matches'])
        group=tuple(row[k] for k in ['guide','policy','mask','screen','feature'])
        for metric in METRICS:
            value=None if row[metric]=='' else float(row[metric]);assert value is None or np.isfinite(value)
            groups[group+(metric,)].append((row['scenario_id'],value))
    assert source_keys==set(itertools.product(guides,policies,scenario_ids,masks,screens,features))
    summaries=[];maxima={}
    for key,values in sorted(groups.items()):
        assert len(values)==54 and len({s for s,x in values})==54
        nums=[x for s,x in values if x is not None];absnums=[abs(x) for x in nums]
        record=dict(zip(['guide','policy','mask','screen','feature','metric'],key),scenarios=54,estimable_scenarios=len(nums),unavailable_scenarios=54-len(nums),minimum_signed=None if not nums else min(nums),maximum_signed=None if not nums else max(nums),median_absolute=None if not nums else float(np.median(absnums)),maximum_absolute=None if not nums else max(absnums))
        summaries.append(record);maxima[key]=record['maximum_absolute']
    assert len(summaries)==4608
    table=Path('docs/tables/full_screened_balance_scenario_summary_20261002.tsv');covtable=Path('docs/tables/full_screened_matching_coverage_20261002.tsv')
    assert not table.exists() and not covtable.exists()
    with table.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(summaries[0]),delimiter='\t');w.writeheader();w.writerows(summaries)
    covtable.write_bytes(paths[1].read_bytes());matrix_rows=list(itertools.product(guides,policies));matrices=[]
    for metric in METRICS:
        matrix=np.asarray([[np.nan if maxima[g,p,'both','n50_c70',feature,metric] is None else maxima[g,p,'both','n50_c70',feature,metric] for feature in features] for g,p in matrix_rows]);matrices.append(matrix)
    maximum=float(np.nanmax(matrices));fig,axes=plt.subplots(2,2,figsize=(15,10),layout='constrained')
    titles=['Matched target–control balance','Retained targets vs all target records','Retained targets vs metadata-matched targets','Retained targets vs quality-eligible targets']
    feature_labels=['Sequence\ndistance','Log positive\nseq. distance','Mean log\nlength','Log length\nasymmetry','Mean\npLDDT','Minimum\npLDDT','Mean low\nconfidence\nfraction','Max. low\nconfidence\nfraction']
    row_labels=[g+' / '+p.replace('alignment_','align. ').replace('envelope_','env. ') for g,p in matrix_rows]
    for ax,metric,title,matrix in zip(axes.ravel(),METRICS,titles,matrices):
        im=ax.imshow(matrix,vmin=0,vmax=maximum,cmap='viridis',aspect='auto');ax.set_title(title,fontsize=12);ax.set_xticks(range(8),feature_labels,fontsize=7);ax.set_yticks(range(8),row_labels,fontsize=8)
        for i,j in itertools.product(range(8),range(8)):
            value=matrix[i,j];ax.text(j,i,'NA' if not np.isfinite(value) else f'{value:.2f}',ha='center',va='center',fontsize=8,color='white' if not np.isfinite(value) or value<maximum*.55 else '#222')
    fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.7,label='Maximum absolute standardized difference across all 54 fixed scenarios')
    fig.suptitle('Post-screen matching balance and target-selection shifts\nBoth masks; ≥50 residues and ≥70% original coverage. Each panel uses its stated SD denominator.',fontsize=14)
    outputs=[]
    for ext in ['png','svg','pdf']:
        target=Path('docs/figures/full_screened_balance_20261002.'+ext);assert not target.exists();fig.savefig(target,dpi=180);outputs.append(target)
    plt.close(fig);ns='{http://www.w3.org/2000/svg}';tree=ET.parse(outputs[1]);text=' '.join(''.join(e.itertext()) for e in tree.iter(ns+'text'))
    for matrix in matrices:
        for value in matrix.ravel():assert ('NA' if not np.isfinite(value) else f'{value:.2f}') in text
    data=dict(metrics=METRICS,rows=matrix_rows,features=features,mask='both',screen='n50_c70',scenario_count=54,maximum_absolute_matrices=[[[None if not np.isfinite(x) else float(x) for x in row] for row in matrix] for matrix in matrices])
    item=ET.Element(ns+'metadata',id='full-screened-balance-values');item.text=json.dumps(data);tree.getroot().append(item);tree.write(outputs[1],encoding='utf-8',xml_declaration=True)
    assert json.loads(ET.parse(outputs[1]).find(ns+'metadata[@id="full-screened-balance-values"]').text)==json.loads(json.dumps(data))
    result=dict(status='published_full_verified_post_screen_balance_and_selection_shift_summary',source_completion=str(a.completion),source_completion_sha256=sha(a.completion),full_source_archive=str(archive),full_source_archive_sha256=sha(archive),original_completion_journals=2,full_balance_rows=62208,full_coverage_rows=7776,summary_rows=4608,plotted_summary_cells=256,plotted_scenarios_per_cell=54,svg_text_and_metadata_values_checked=True,
        main_screen_maximum_absolute_by_metric={metric:float(np.nanmax(matrix)) for metric,matrix in zip(METRICS,matrices)},source_hashes={str(p):sha(p) for p in [a.completion,archive,rp,ap,*paths,Path(__file__)]},artifacts={str(p):sha(p) for p in [table,covtable,*outputs]},scientific_eligibility=False,scope='Full62208balance rows/7776coverage strata retained with4608all-mask/screen/guide/policy/feature/metric summaries acrossall54scenarios. Report signedextrema/absolute medians/maxima/estimable/missingcounts; figureshowsall256bothmaskn50c70 maxima, never favorable matching scenario. Magnitudes use distinct source denominators permetric; no aggregation asindependent events orweightedbiologicaleffects. Sourcefullhash/journal refresh and actualfigurevisualcheck separate. Sourcebalance doesnot proveunbiased sampling or calibratedcausal/duplicationeffects;all8aims incomplete;GPUpaused.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}),flush=True)


if __name__=='__main__':main()
