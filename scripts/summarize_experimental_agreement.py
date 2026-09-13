#!/usr/bin/env python3
"""Summarize descriptive agreement by protein and fixed chain/model cohorts."""
import argparse,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--comparisons',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use an immutable output')
    r=checked_receipt(a.comparisons);rows=read_table(a.comparisons/'comparisons.tsv');excluded=read_table(a.comparisons/'exclusions.tsv')
    if len(rows)!=r['accepted_rows'] or len(excluded)!=r['excluded_rows']:raise ValueError('Receipt count mismatch')
    key=lambda x:tuple(x[k] for k in ['entry_id','entity_id','deposited_model','label_asym_id','predicted_model_id'])
    seen=set();bycut=defaultdict(dict);groups=defaultdict(list)
    for accepted,items in [(True,rows),(False,excluded)]:
        for row in items:
            cut=int(row['predicted_plddt_cutoff']);identity=key(row)
            if (identity,cut) in seen or cut not in (0,70,90):raise ValueError('Duplicate/invalid comparison')
            seen.add((identity,cut));n=int(row['matched_residues']);length=int(row['sequence_length'])
            if accepted!=(n>=50 and 2*n>=length):raise ValueError('Eligibility inconsistent')
            if accepted:bycut[cut][identity]=row;groups[cut,row['predicted_model_id']].append(row)
    if any({cut for k,cut in seen if k==identity}!={0,70,90} for identity,_ in seen):raise ValueError('Incomplete threshold grid')
    protein=[]
    for (cut,model),items in sorted(groups.items()):
        protein.append({'predicted_model_id':model,'predicted_plddt_cutoff':cut,'chain_model_comparisons':len(items),'entries':len({r['entry_id'] for r in items}),'median_CA_RMSD_angstrom':float(np.median([float(r['ca_superposition_rmsd_angstrom']) for r in items])),'experimental_methods':';'.join(sorted({r['experimental_method'] for r in items}))})
    cohorts=[]
    for cut in [70,90]:
        ids=set(bycut[cut]);baseline=bycut[0]
        if not ids<=baseline.keys():raise ValueError('Higher threshold lacks baseline')
        cohorts.append({'cutoff':cut,'chain_model_comparisons':len(ids),'distinct_predicted_proteins':len({i[-1] for i in ids}),'baseline_median_CA_RMSD_angstrom':float(np.median([float(baseline[i]['ca_superposition_rmsd_angstrom']) for i in ids])),'filtered_median_CA_RMSD_angstrom':float(np.median([float(bycut[cut][i]['ca_superposition_rmsd_angstrom']) for i in ids])),'interpretation':'Identical chain/model cohort; residue subsets differ by predicted focal confidence. No causal accuracy improvement or independent samples implied.'})
    a.output.mkdir(parents=True);write_table(a.output/'protein_summary.tsv',protein);write_table(a.output/'matched_cohort_summary.tsv',cohorts)
    result={'status':'complete_descriptive_protein_and_matched_cohort_summary','comparison_receipt_sha256':sha(a.comparisons/'receipt.json'),'script_sha256':sha(Path(__file__)),'unique_chain_model_targets':len(seen)//3,'threshold_grid_rows':len(seen),'equal_protein_median_of_within_protein_median_RMSD':{str(cut):float(np.median([r['median_CA_RMSD_angstrom'] for r in protein if r['predicted_plddt_cutoff']==cut])) for cut in [0,70,90]},'matched_cohorts':cohorts,'interpretation':'Descriptive summaries of the retrieval-order partial snapshot. Equal-protein medians prevent prolific depositions from receiving extra between-protein weight but do not make proteins representative or independent. Each protein summary still weights all its chain/model observations equally. Experimental quality, context and training independence remain unresolved.','artifacts':{f.name:sha(f) for f in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
