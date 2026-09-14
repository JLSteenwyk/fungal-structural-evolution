#!/usr/bin/env python3
"""Match all baseline/FCS conditional coefficients without testing their difference."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['baseline','sensitivity','output']:
        p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    receipts={k:checked_receipt(getattr(a,k)) for k in ['baseline','sensitivity']}
    frames={k:pd.read_csv(getattr(a,k)/'focal_coefficients.tsv',sep='\t',float_precision='round_trip') for k in receipts}
    keys=['model_id','alphabet_model','rate_model','rsa_scale','controls','term']
    for k,r in receipts.items():
        if r['status']!='complete_exploratory_conditional_site_coupling_fits' or r['fits']!=24 or r['focal_tests']!=72 or frames[k].duplicated(keys).any():
            raise ValueError('Incomplete or duplicate planned fit grid')
    if any(receipts['baseline'][k]!=receipts['sensitivity'][k] for k in ['sites','markers']):
        raise ValueError('Site/marker counts differ')
    joined=frames['baseline'].merge(frames['sensitivity'],on=keys,how='outer',suffixes=('_baseline','_sensitivity'),validate='one_to_one',indicator=True)
    if len(joined)!=72 or not joined._merge.eq('both').all():raise ValueError('Model identity differs')
    joined=joined.drop(columns='_merge')
    joined['coefficient_change']=joined.coefficient_sensitivity-joined.coefficient_baseline
    joined['coefficient_sign_changed']=np.sign(joined.coefficient_sensitivity)!=np.sign(joined.coefficient_baseline)
    summaries=[]
    for term,g in joined.groupby('term',sort=True):
        summaries.append({'term':term,'models':len(g),'sign_changes':int(g.coefficient_sign_changed.sum()),
                          'baseline_bh_q_below_0_05':int((g.bh_q_across_all_72_focal_tests_baseline<.05).sum()),
                          'sensitivity_bh_q_below_0_05':int((g.bh_q_across_all_72_focal_tests_sensitivity<.05).sum()),
                          'sensitivity_coefficient_minimum':float(g.coefficient_sensitivity.min()),
                          'sensitivity_coefficient_maximum':float(g.coefficient_sensitivity.max()),
                          'maximum_absolute_coefficient_change':float(abs(g.coefficient_change).max()),
                          'sensitivity_intervals_containing_zero':int(((g.ci95_lower_sensitivity<=0)&(g.ci95_upper_sensitivity>=0)).sum())})
    a.output.mkdir(parents=True)
    joined.to_csv(a.output/'matched_coefficients.tsv',sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(a.output/'term_summary.tsv',sep='\t',index=False)
    result={'status':'complete_matched_fcs_conditional_coupling_comparison','coefficient_pairs':72,
            'source_receipts':{k:sha(getattr(a,k)/'receipt.json') for k in receipts},
            'script_sha256':sha(Path(__file__)),'summaries':summaries,
            'artifacts':{x.name:sha(x) for x in a.output.iterdir()},
            'interpretation':'Descriptive differences between matched baseline and omission fits using the same specifications. Each fit set uses its own original 72-test BH family; threshold crossing is not a test of coefficient change. Models, sites and datasets are dependent. Conditional marker-cluster intervals do not propagate rate/tree/prediction uncertainty. FCS omission is a sensitivity analysis, not proof of contamination or robust lineage-wide coupling.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
