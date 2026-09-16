#!/usr/bin/env python3
"""Summarize paired predictor controls with hierarchical medians and equal protein weight."""
import argparse
import json
from collections import defaultdict
from pathlib import Path
from statistics import median
from audit_joint_path_uncertainty import checked, rows, sha, quantile
from prepare_paired_phylogenetic_inputs import write_table

METRICS=[label+'_'+metric for label in ['af_experiment','esm_experiment','af_esm'] for metric in ['rmsd_angstrom','local_mean_absolute_difference_angstrom','local_rms_difference_angstrom']]
DELTAS=['paired_delta_'+metric for metric in ['rmsd_angstrom','local_mean_absolute_difference_angstrom','local_rms_difference_angstrom']]


def aggregate(data,keys):
    groups=defaultdict(list)
    for row in data:groups[tuple(row[k] for k in keys)].append(row)
    output=[]
    for identity,items in sorted(groups.items()):
        row=dict(zip(keys,identity));row['contributing_units']=len(items)
        for metric in METRICS+DELTAS+['fraction_full_sequence']:
            values=[float(x[metric]) for x in items if x[metric]!='']
            row[metric]=median(values) if values else ''
        output.append(row)
    return output


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['comparisons','audit','crosswalk','output']:ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=checked(a.comparisons);c=checked(a.crosswalk);audit=json.loads(a.audit.read_text())
    if audit['status']!='passed_full_experimental_control_geometry_readback' or audit['comparison_receipt_sha256']!=sha(a.comparisons/'receipt.json') or r['source_receipts']['crosswalk']!=sha(a.crosswalk/'receipt.json'):
        raise ValueError('Matching full geometry audit required')
    data=list(rows(a.comparisons/'comparisons.tsv'));excluded=list(rows(a.comparisons/'exclusions.tsv'))
    if len(data)!=r['accepted_rows'] or len(excluded)!=r['excluded_rows']:raise ValueError('Row count mismatch')
    for row in data:
        for metric in ['rmsd_angstrom','local_mean_absolute_difference_angstrom','local_rms_difference_angstrom']:
            af,esm=row['af_experiment_'+metric],row['esm_experiment_'+metric]
            if (af=='')!=(esm==''):raise ValueError('Unpaired metric availability')
            row['paired_delta_'+metric]=float(esm)-float(af) if af!='' else ''
    keys=['sequence_sha256','joint_predicted_plddt_cutoff']
    models=aggregate(data,keys+['entry_id','deposited_model'])
    entries=aggregate(models,keys+['entry_id']);proteins=aggregate(entries,keys)
    universe={x['sequence_sha256'] for x in rows(a.crosswalk/'crosswalk.tsv')}
    bykey={(x['sequence_sha256'],x['joint_predicted_plddt_cutoff']):x for x in proteins}
    if len(bykey)!=len(proteins) or len(universe)!=c['control_sequences']:raise ValueError('Protein universe differs')
    dispositions=[]
    for sequence in sorted(universe):
        for cutoff in ['0','70','90']:
            accepted=[x for x in data if x['sequence_sha256']==sequence and x['joint_predicted_plddt_cutoff']==cutoff]
            rejected=[x for x in excluded if x['sequence_sha256']==sequence and x['joint_predicted_plddt_cutoff']==cutoff]
            dispositions.append({'sequence_sha256':sequence,'joint_predicted_plddt_cutoff':cutoff,'accepted_chain_model_rows':len(accepted),'excluded_chain_model_rows':len(rejected),'accepted_entries':len({x['entry_id'] for x in accepted}),'status':'included' if accepted else 'no_eligible_chain_model'})
    common=set.intersection(*[{x['sequence_sha256'] for x in proteins if x['joint_predicted_plddt_cutoff']==cutoff} for cutoff in ['0','70','90']])
    summaries=[]
    for cohort in ['all_available','same_proteins_all_thresholds']:
        for cutoff in ['0','70','90']:
            subset=[x for x in proteins if x['joint_predicted_plddt_cutoff']==cutoff and (cohort=='all_available' or x['sequence_sha256'] in common)]
            for metric in METRICS+DELTAS:
                values=[float(x[metric]) for x in subset if x[metric]!='']
                summaries.append({'cohort':cohort,'joint_predicted_plddt_cutoff':cutoff,'metric':metric,'proteins':len(subset),'proteins_with_metric':len(values),'median':median(values) if values else '', 'q25':quantile(values,.25) if values else '', 'q75':quantile(values,.75) if values else '', 'positive':sum(x>1e-8 for x in values) if metric in DELTAS else '', 'negative':sum(x< -1e-8 for x in values) if metric in DELTAS else '', 'approximately_zero':sum(abs(x)<=1e-8 for x in values) if metric in DELTAS else ''})
    a.output.mkdir(parents=True)
    for name,values in [('deposited_model_summary.tsv',models),('entry_summary.tsv',entries),('protein_summary.tsv',proteins),('protein_disposition.tsv',dispositions),('cohort_summary.tsv',summaries)]:write_table(a.output/name,values)
    result={'status':'complete_hierarchical_experimental_predictor_control_summary','control_sequences':len(universe),'protein_threshold_rows':len(proteins),'disposition_rows':len(dispositions),'same_proteins_all_thresholds':len(common),'comparison_receipt_sha256':sha(a.comparisons/'receipt.json'),'audit_sha256':sha(a.audit),'script_sha256':sha(Path(__file__)),'interpretation':'Median chains within deposited model, median models within entry, median entries within canonical protein; then equal-protein descriptive medians/quartiles. Paired deltas are ESM-experiment minus AF-experiment at the chain level before hierarchical aggregation; positive favors lower AF discrepancy. Matched protein cohorts still use threshold-dependent residues and may use different entries. No independent-protein assumption, significance test, calibrated accuracy, training-independence or causal claim.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
