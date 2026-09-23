#!/usr/bin/env python3
"""Join all codon-case diagnostics without declaring selection eligibility."""
from collections import Counter
import csv
import json
from pathlib import Path
from catalog_whole_proteome_structures import sha

SOURCES={
 'information':'metadata/genus_codon_tree_information.tsv',
 'trees':'metadata/genus_codon_tree_full_case_review.tsv',
 'fits':'metadata/genus_mg94_full_case_audit.tsv',
 'profiles':'metadata/genus_mg94_branch_parameter_profile_audit_case_summary.tsv',
 'fcs':'metadata/fcs_codon_omission_disposition.tsv',
 'restarts':'metadata/genus_mg94_flagged_restart_case_summary.tsv'}


def main():
    hashes={p:sha(p) for p in SOURCES.values()};tables={}
    for name,path in SOURCES.items():
        with Path(path).open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        tables[name]={r['case_id']:r for r in rows}
        if len(rows)!=len(tables[name]):raise ValueError('Duplicate case')
    info=tables['information'];fitted=set(tables['fits'])
    expected={k for k,r in info.items() if r['status']=='ready_for_supported_tree_diagnostic'}
    if fitted!=expected or any(set(tables[k])!=fitted for k in ['trees','profiles','fcs']):raise ValueError('Inconsistent full diagnostic case grid')
    flagged={k for k,r in tables['profiles'].items() if float(r['best_grid_exceeds_unconstrained_by'])>1e-5}
    if set(tables['restarts'])!=flagged:raise ValueError('Restart scope differs from profile flags')
    rows=[];counts=Counter()
    for case,r in sorted(info.items()):
        reasons=[];has_fit=case in fitted
        if not has_fit:reasons.append('insufficient_sequence_information')
        if r['marker_copy_caveat']!='none_recorded':reasons.append('copy_reconciliation_required')
        if has_fit:
            t=tables['trees'][case];fit=tables['fits'][case];prof=tables['profiles'][case];fcs=tables['fcs'][case]
            if any(x['marker_copy_caveat']!=r['marker_copy_caveat'] for x in [t,fit,prof]):raise ValueError('Copy status disagrees')
            if not int(r['taxa'])==int(t['taxa'])==int(fit['taxa'])==int(fcs['original_taxa']):raise ValueError('Taxon counts disagree')
            if int(r['nucleotide_columns'])!=3*int(fit['codons']) or int(t['nucleotide_columns'])!=int(r['nucleotide_columns']):raise ValueError('Alignment dimensions disagree')
            if fcs['status']=='below_existing_four_taxon_minimum':reasons.append('fcs_omission_below_four_taxa')
            elif fcs['status']!='unchanged_by_declared_fcs_omission':raise ValueError('Unknown FCS status')
            for field,label in [('saturated_pairwise_distance_warning_lines','tree_pairwise_saturation_warning'),('parameter_boundary_warning_lines','tree_parameter_boundary_warning'),('nni_convergence_warning_lines','tree_nni_convergence_warning'),('other_warning_lines','other_tree_warning')]:
                if int(t[field]):reasons.append(label)
            if int(t['long_edges_ge10']):reasons.append('very_long_nucleotide_tree_edge')
            if int(t['near_zero_edges_le1e_5']):reasons.append('near_zero_nucleotide_tree_edge')
            if fit['numerical_review_flags']:reasons.append('saved_fit_numerical_review')
            if case in flagged:reasons.append('multistart_refit_identifiability_unresolved')
        else:t=fit=prof=fcs={}
        for reason in reasons:counts[reason]+=1
        rows.append(dict(case_id=case,taxa=r['taxa'],nucleotide_columns=r['nucleotide_columns'],distinct_aligned_sequences=r['distinct_aligned_sequences'],informative_nucleotide_columns=r['parsimony_informative_nucleotide_columns'],diagnostic_fit_completed=has_fit,copy_caveat=r['marker_copy_caveat'],fcs_omission_status=fcs.get('status','not_fitted'),remaining_taxa_after_fcs_omission=fcs.get('remaining_taxa',''),target_original_dS_review_label=prof.get('original_dS_review_label',''),multistart_refit_completed=case in flagged,case_specific_review_flags=';'.join(reasons),selection_eligibility='not_established',remaining_general_requirements='saturation_and_identifiability;alignment_and_copy_adequacy;selection_model_and_test_design'))
    out=Path('results/cds/codon-analysis-readiness-20260923-v1')
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);table=out/'case_readiness.tsv'
    with table.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    for path,h in hashes.items():
        if sha(path)!=h:raise ValueError('Changed diagnostic input')
    result={'status':'complete_full_codon_case_review_ledger','input_cases':len(info),'fitted_cases':len(fitted),'unfitted_cases':len(info)-len(fitted),'multistart_review_cases':len(flagged),'case_specific_flag_counts':dict(counts),'cases_with_case_specific_flags':sum(bool(r['case_specific_review_flags']) for r in rows),'cases_without_listed_case_specific_flags':sum(not r['case_specific_review_flags'] for r in rows),'selection_eligible_cases':0,'sources':hashes,'script_sha256':sha(__file__),'artifacts':{'case_readiness.tsv':sha(table)},'scope':'Exact joins of completed diagnostic summaries, with case grids, copy flags, taxon counts and alignment dimensions checked. Flags overlap and are review triggers, not validated exclusion thresholds. Absence of listed flags does not establish eligibility. No new selection test, saturation certificate or biological rate conclusion.'}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
