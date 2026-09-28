"""Check real scalar records plus missing/tampered/failed dispositions."""
import copy
import json
from pathlib import Path
import pandas as pd
from ancestral_chain_attempt import sha,write_json
from audit_selected_matched_intervals import check_fit,NAMES

fits=json.loads(Path('results/model_validation/selected-matched-interval-checks-20260928-v1/real_input_results.json').read_text())
frame=pd.read_parquet('results/structural_comparisons/refined-working-model-grid-export-20260928-v2/unique_fits.parquet').set_index(['fit_input_id','tree'],drop=False)
for fit in fits:
    row=frame.loc[(fit['fit_input_id'],fit['tree'])].to_dict()
    assert len(check_fit(fit,row))==5
for change in ['endpoint','estimate','missing','source']:
    bad=copy.deepcopy(fit)
    if change=='endpoint':bad['intervals']['intercept']['upper']+=1
    if change=='estimate':bad['intervals']['intercept']['estimate']+=1
    if change=='missing':del bad['intervals']['intercept']
    if change=='source':bad['selected_source_sha256']='changed'
    try:check_fit(bad,row)
    except AssertionError:pass
    else:raise AssertionError('Invalid record accepted: '+change)
failed={k:fit[k] for k in ['fit_input_id','tree','selection','selected_source_sha256']}
failed.update(status='interval_computation_requires_review',intervals={},error_type='ValueError',error='fixture')
assert len(check_fit(failed,row))==5
omitted=copy.deepcopy(fit);omitted_row=copy.deepcopy(row)
name=NAMES[-1];omitted_row['selected_coefficient_'+name]=None
omitted['intervals'][name]=dict(status='omitted_constant_covariate',estimate=None)
assert len(check_fit(omitted,omitted_row))==5
result=dict(status='passed_real_intervals_and_audit_rejection_fixtures',real_fits=len(fits),
    rejected=['changed_endpoint','changed_estimate','missing_coefficient','changed_source'],
    failed_fit_retains_five_dispositions=True,omission_retained=True,
    pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/audit_selected_matched_intervals.py')]})
write_json(Path('metadata/selected_matched_interval_audit_checks_20260928.json'),result)
print(json.dumps(result))
