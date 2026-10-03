#!/usr/bin/env python3
"""Export complete closed matched-source branch controls and descriptive figures."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from ancestral_chain_attempt import sha
from background_measurement_union_sources import closed_source
from matched_predictor_branch_fits import load
from matched_predictor_branch_inputs import verify

MODELS=['af','af_empirical','llm']
LABELS=['AF + G4','AF + empirical frequencies + G4','LLM + G4']
REVIEW_THRESHOLD=1e-6


def table(path,records):
    with path.open('x') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(records)
    with path.open() as f:
        actual=list(csv.DictReader(f,delimiter='\t'))
    assert actual==[{k:'' if v is None else str(v) for k,v in r.items()} for r in records]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists() and not a.receipt.exists()
    plan_path=Path('metadata/matched_predictor_branch_fits_plan_20261003_v1.json');plan=json.loads(plan_path.read_text())
    source,bindings=load(plan,plan_path)
    completion=closed_source(plan['completion'],'complete_verified_full_matched_predictor_native_point_fits',
        'complete_verified_full_matched_predictor_native_point_fits_archive',2,bindings)
    assert completion['native_roles']==931 and completion['unresolved_inputs']==0
    bindings[str(Path(__file__))]=sha(__file__);root=Path(plan['output']);a.output.mkdir(parents=True)
    pairs=[];summaries=[]
    for key,config in sorted(source['configs'].items()):
        def points(role):
            record=json.loads((root/'roles'/(key+'-'+role+'.json')).read_text())
            assert record['status']=='native_point_output_integrity_checked_not_model_qualified'
            return {r['split_mask_hex']:r for r in record['point_estimate']['branches']}
        aa=points('aa')
        for model in MODELS:
            af,esm=points('AlphaFold_'+model),points('ESMFold_'+model);assert set(aa)==set(af)==set(esm)
            for split in sorted(aa,key=lambda x:int(x,16)):
                baseline,first,second=aa[split]['length'],af[split]['length'],esm[split]['length']
                assert aa[split]['internal']==af[split]['internal']==esm[split]['internal']
                difference=second-first;total=second+first;normalized=difference/total if total>0 else None
                assert normalized is None or -1<=normalized<=1
                pairs.append(dict(input_id=key,marker=config['marker'],model=model,split_mask_hex=split,
                    branch_class='projected_internal' if aa[split]['internal'] else 'terminal',taxa=len(config['taxa']),columns=len(config['columns']),
                    aa_length=baseline,alphafold_length=first,esmfold_length=second,esm_minus_af=difference,
                    symmetric_source_difference=normalized,aa_at_or_below_review_threshold=baseline<=REVIEW_THRESHOLD,
                    af_at_or_below_review_threshold=first<=REVIEW_THRESHOLD,esm_at_or_below_review_threshold=second<=REVIEW_THRESHOLD,
                    branch_unit='expected_state_substitutions_per_site',scientific_eligibility=False))
    assert len(pairs)==6585 and sum(r['branch_class']=='projected_internal' for r in pairs)==2694
    table(a.output/'all_matched_branch_points.tsv',pairs)
    for model in MODELS:
        for kind in ['projected_internal','terminal']:
            selected=[r for r in pairs if r['model']==model and r['branch_class']==kind]
            normalized=[abs(r['symmetric_source_difference']) for r in selected if r['symmetric_source_difference'] is not None]
            summaries.append(dict(model=model,branch_class=kind,conditional_branch_pairs=len(selected),
                markers=len({r['marker'] for r in selected}),inputs=len({r['input_id'] for r in selected}),
                median_af_length=float(np.median([r['alphafold_length'] for r in selected])),
                median_esm_length=float(np.median([r['esmfold_length'] for r in selected])),
                median_absolute_source_length_difference=float(np.median([abs(r['esm_minus_af']) for r in selected])),
                median_absolute_symmetric_source_difference=float(np.median(normalized)),
                q90_absolute_symmetric_source_difference=float(np.quantile(normalized,.9)),
                af_at_or_below_review_threshold=sum(r['af_at_or_below_review_threshold'] for r in selected),
                esm_at_or_below_review_threshold=sum(r['esm_at_or_below_review_threshold'] for r in selected),
                both_source_lengths_zero=sum(r['symmetric_source_difference'] is None for r in selected),
                review_threshold=REVIEW_THRESHOLD,scientific_eligibility=False))
    table(a.output/'model_source_sensitivity_summary.tsv',summaries)
    with gzip.open(source['root']/'comparison_cases.jsonl.gz','rt') as f:cases=[json.loads(line) for line in f]
    coverage=[]
    for vi,view in enumerate(source['axes']['views']):
        selected=[r for r in cases if r['view_index']==vi];assert len(selected)==125
        counts=Counter(r['status'] for r in selected)
        coverage.append(dict(view_index=vi,cohort=view['cohort'],view=view['view'],marker_slots=125,
            ready_markers=counts['ready_for_matched_predictor_branch_fits'],insufficient_markers=counts['insufficient_joint_observation_taxa'],
            original_internal_branches=len(view['branches']),original_internal_marker_branch_slots=sum(r['original_internal_branches'] for r in selected),
            unique_original_branch_projections=sum(r['branch_status_counts']['unique_internal_projection'] for r in selected),
            shared_original_branch_projections=sum(r['branch_status_counts']['shared_internal_projection'] for r in selected)))
    table(a.output/'all_view_control_coverage.tsv',coverage)
    marker_counts=[]
    for marker in source['axes']['markers']:
        selected=[r for r in cases if r['marker']==marker];assert len(selected)==70
        native=[r for r in pairs if r['marker']==marker]
        marker_counts.append(dict(marker=marker,tree_views=70,ready_views=sum(r['input_id'] is not None for r in selected),
            insufficient_views=sum(r['input_id'] is None for r in selected),unique_native_inputs=len({r['input_id'] for r in native}),
            matched_taxa_union=len({t for r in selected for t in r['taxa']}),conditional_internal_branch_model_pairs=sum(r['branch_class']=='projected_internal' for r in native)))
    table(a.output/'all_marker_control_coverage.tsv',marker_counts)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(2,3,figsize=(12.5,8),constrained_layout=True)
    for col,(model,label) in enumerate(zip(MODELS,LABELS)):
        selected=[r for r in pairs if r['model']==model and r['branch_class']=='projected_internal']
        first=np.asarray([r['alphafold_length'] for r in selected]);second=np.asarray([r['esmfold_length'] for r in selected])
        valid=(first>0)&(second>0);upper=max(first.max(),second.max())*1.2
        lower=min(first[valid].min(),second[valid].min())*.7
        ax=axes[0,col];ax.scatter(first[valid],second[valid],s=12,alpha=.38,color='#126782',linewidths=0,rasterized=True)
        ax.plot([lower,upper],[lower,upper],color='.5',lw=1,ls='--');ax.set_xscale('log');ax.set_yscale('log')
        ax.set_xlim(lower,upper);ax.set_ylim(lower,upper);ax.set_title(label);ax.set_xlabel('AlphaFold structural branch length')
        if col==0:ax.set_ylabel('ESMFold structural branch length')
        ax.text(.04,.96,f'{len(selected)} conditional internal branch pairs',transform=ax.transAxes,va='top',fontsize=9)
        ax=axes[1,col]
        x=np.asarray([r['aa_length'] for r in selected]);y=np.asarray([np.nan if r['symmetric_source_difference'] is None else r['symmetric_source_difference'] for r in selected])
        ok=(x>0)&np.isfinite(y);ax.scatter(x[ok],y[ok],s=12,alpha=.38,color='#bb580d',linewidths=0,rasterized=True)
        ax.axhline(0,color='.5',lw=1);ax.set_xscale('log');ax.set_ylim(-1.05,1.05)
        ax.set_xlabel('Amino-acid branch length (LG + F + G4)')
        if col==0:ax.set_ylabel('(ESMFold − AlphaFold) / sum')
    fig.suptitle('Predictor sensitivity with identical proteins, observations and fixed topologies',fontsize=14)
    for ext in ['png','pdf']:fig.savefig(a.output/('matched_predictor_branch_sensitivity.'+ext),dpi=180)
    plt.close(fig);verify(bindings)
    result=dict(status='complete_descriptive_closed_matched_predictor_branch_summary',source_hashes=bindings,
        compared_native_branch_pairs=len(pairs),conditional_internal_branch_pairs=2694,models=MODELS,markers=71,taxon_entries=21,
        full_scope_views=70,full_scope_marker_slots=125,view_coverage_rows=70,marker_coverage_rows=125,
        model_summaries=summaries,artifacts={str(p):sha(p) for p in a.output.iterdir()},
        scientific_eligibility=False,scope='Every closed native branch point joined across identical AA/predictor inputs. All6585conditional paired branch values including2694internal pairs retained, plus all70views/125marker slots. Exact TSV readback, full upstream hash proofs and original closures verified. Figures/descriptive unweighted quantiles condition on133input/topology combinations and reuse71markers/21fungal taxa; these are correlated alternatives, not independent replicates or unbiased lineage samples. Raw1e-6 review threshold flags tiny native branches without claiming IQTree default/boundary or replacing any value. No significance, calibrated uncertainty, accepted evolutionary acceleration, physical displacement, ecological/duplication effect or biological aim completion.')
    with a.receipt.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']},indent=2),flush=True)


if __name__=='__main__':main()
