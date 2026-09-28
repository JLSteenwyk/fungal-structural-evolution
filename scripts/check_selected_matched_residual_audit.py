"""Check nested replay comparison rejects changed diagnostics and lost statuses."""
import copy
from pathlib import Path
from ancestral_chain_attempt import sha,write_json
from audit_selected_matched_residuals import compare
example=dict(status='descriptive_marginal_residual_diagnostics',summaries=dict(mean=.1,quantiles=[-.4,.2,1.3]),covariate_bins=[dict(records=12,mean_standardized_residual=.2)])
compare(copy.deepcopy(example),example)
for field in ['mean','quantile','bin','status']:
    bad=copy.deepcopy(example)
    if field=='mean':bad['summaries']['mean']+=.01
    elif field=='quantile':bad['summaries']['quantiles'][0]+=.01
    elif field=='bin':bad['covariate_bins']=[]
    else:bad['status']='changed'
    try:compare(bad,example)
    except AssertionError:pass
    else:raise AssertionError('Changed '+field+' accepted')
failed=dict(status='residual_computation_requires_review',error_type='ValueError',error='fixture')
compare(failed,failed.copy())
write_json(Path('metadata/selected_matched_residual_audit_checks_20260928.json'),dict(status='passed_nested_diagnostic_rejection_and_failure_preservation_checks',rejected_changes=4,failed_status_preserved=True,pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/audit_selected_matched_residuals.py')]},scope='Replay comparison checks; evaluator real-data and dense projection checks are recorded separately.'))
print('Passed nested residual audit checks')
