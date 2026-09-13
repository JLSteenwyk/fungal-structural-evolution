#!/usr/bin/env python3
"""Summarize descriptive comparisons without treating repeated pairs as independent."""
import argparse,json
from pathlib import Path
import pandas as pd
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--comparisons',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new summary directory')
    checked_receipt(a.comparisons)
    df=pd.read_csv(a.comparisons/'pairwise_metrics.tsv',sep='\t');ex=pd.read_csv(a.comparisons/'pairwise_exclusions.tsv',sep='\t');rows=[]
    for cutoff,sub in df.groupby('plddt_cutoff'):
        rows.append({'plddt_cutoff':int(cutoff),'accepted_pairs':len(sub),'excluded_pairs':int((ex.plddt_cutoff==cutoff).sum()),'markers':int(sub.marker.nunique()),'median_compared_residues':float(sub.compared_residues.median()),'median_uncorrected_sequence_difference':float(sub.uncorrected_sequence_difference.median()),'median_ca_rmsd_angstrom':float(sub.ca_superposition_rmsd_angstrom.median()),'median_local_distance_change_angstrom':float(sub.local_distance_mean_absolute_change_angstrom.median())})
    a.output.mkdir(parents=True);table=a.output/'threshold_summary.tsv';pd.DataFrame(rows).to_csv(table,sep='\t',index=False)
    result={'status':'complete_descriptive_geometry_summary','source_receipt_sha256':sha(a.comparisons/'receipt.json'),'script_sha256':sha(Path(__file__)),'pandas_version':pd.__version__,'artifacts':{table.name:sha(table)},'interpretation':'Pair-weighted descriptive medians. Pairs share taxa and families; threshold strata differ in accepted pairs and residues. No independent-observation confidence interval, rate, predictor accuracy or causal threshold effect is inferred.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
