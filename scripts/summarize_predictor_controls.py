#!/usr/bin/env python3
"""Summarize predictor differences with explicit changing and shared cohorts."""
import argparse
import json
import math
import statistics
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def median(rows,key):
    values=[float(r[key]) for r in rows if r[key]!='']
    return statistics.median(values) if values else ''


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--comparisons',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new summary output')
    receipt=checked_receipt(a.comparisons);rows=read_table(a.comparisons/'comparisons.tsv');lookup={}
    for r in rows:
        key=(r['sequence_id'],int(r['plddt_cutoff']))
        if key in lookup:raise ValueError('Duplicate control comparison')
        lookup[key]=r
        if r['status']=='compared':
            for name in ['ca_superposition_rmsd_angstrom','local_distance_mean_absolute_change_angstrom','local_distance_rms_change_angstrom']:
                if r[name]!='' and (not math.isfinite(float(r[name])) or float(r[name])<0):raise ValueError('Invalid geometry metric')
            counts=[int(r[f'pae{p}_local_pairs']) for p in [5,10,15]]+[int(r['local_distance_pairs'])]
            if counts!=sorted(counts):raise ValueError('Non-nested PAE confidence sets')
    ids={k[0] for k in lookup}
    if set(lookup)!={(s,c) for s in ids for c in [0,70,90]} or len(ids)!=receipt['unique_sequences']:raise ValueError('Incomplete control grid')
    summaries=[]
    for cutoff in [0,70,90]:
        subset=[r for r in rows if int(r['plddt_cutoff'])==cutoff and r['status']=='compared']
        summaries.append({'plddt_cutoff':cutoff,'total_controls':len(ids),'compared_controls':len(subset),'coverage_excluded_controls':len(ids)-len(subset),
            'median_ca_rmsd_angstrom':median(subset,'ca_superposition_rmsd_angstrom'),'median_local_distance_change_angstrom':median(subset,'local_distance_mean_absolute_change_angstrom'),
            'controls_with_pae10_local_pairs':sum(int(r['pae10_local_pairs'])>0 for r in subset),
            'median_pae10_local_distance_change_angstrom':median(subset,'pae10_local_mean_absolute_change_angstrom')})
    paired=[]
    for low,high in [(0,70),(0,90),(70,90)]:
        shared=sorted(s for s in ids if lookup[(s,low)]['status']==lookup[(s,high)]['status']=='compared')
        before=[lookup[(s,low)] for s in shared];after=[lookup[(s,high)] for s in shared]
        paired.append({'lower_plddt_cutoff':low,'higher_plddt_cutoff':high,'same_controls':len(shared),
            'lower_cutoff_median_rmsd_angstrom':median(before,'ca_superposition_rmsd_angstrom'),'higher_cutoff_median_rmsd_angstrom':median(after,'ca_superposition_rmsd_angstrom'),
            'median_within_control_rmsd_change_angstrom':statistics.median(float(lookup[(s,high)]['ca_superposition_rmsd_angstrom'])-float(lookup[(s,low)]['ca_superposition_rmsd_angstrom']) for s in shared) if shared else ''})
    a.output.mkdir(parents=True);write_table(a.output/'threshold_summary.tsv',summaries);write_table(a.output/'same_control_threshold_changes.tsv',paired)
    r={'status':'complete_predictor_control_descriptive_summary','comparison_receipt_sha256':sha(a.comparisons/'receipt.json'),'script_sha256':sha(Path(__file__)),'unique_controls':len(ids),'comparison_rows':len(rows),
       'interpretation':'Threshold summaries change both protein cohort and compared residues. Shared-cohort summaries hold proteins constant but not residue positions. No accuracy, causality, rate, independence or statistical-significance claim. Source controls intentionally stratify lineage, length and AF confidence; medians are not population-wide estimates.',
       'artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(summaries,indent=2));print(json.dumps(paired,indent=2))

if __name__=='__main__':main()
