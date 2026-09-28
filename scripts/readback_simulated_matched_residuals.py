"""Read back all residual-reference counts and empirical summary arithmetic."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from ancestral_chain_attempt import sha,write_json

root=Path('results/model_validation/matched-simulation-residual-reference-20260928-v1')
r=json.loads((root/'receipt.json').read_text())
assert r['status']=='complete_existing_simulation_residual_summaries_pending_readback'
for name,h in r['artifacts'].items():assert sha(root/name)==h
frame=pd.read_parquet(root/'replicate_diagnostics.parquet')
table=pd.read_csv(root/'diagnostic_reference_summary.tsv',sep='\t')
assert len(frame)==31968 and not frame.duplicated(['case','replicate','mode']).any()
assert len(table)==256 and not table.duplicated(['case','mode','metric']).any()
assert frame.status.value_counts().to_dict()==r['status_counts']
for (case,mode),group in frame.groupby(['case','mode']):
    assert len(group)==999 and set(group.replicate)==set(range(999))
    assert mode in ['generating_covariance','fitted_covariance']
for row in table.to_dict('records'):
    group=frame[(frame.case==row['case'])&(frame['mode']==row['mode'])]
    valid=group[group.status=='descriptive_marginal_residual_diagnostics']
    assert row['attempted']==999 and row['available']==len(valid) and row['unresolved']==999-len(valid)
    values=valid[row['metric']].to_numpy()
    assert np.isfinite(values).all()
    if len(values):
        expected=np.r_[np.mean(values),np.quantile(values,[.025,.5,.975])]
        np.testing.assert_allclose([row[k] for k in ['mean','q025','median','q975']],expected,rtol=1e-12,atol=1e-12)
    else:assert all(pd.isna(row[k]) for k in ['mean','q025','median','q975'])
assert int(frame.qualified_by_continuation.sum())==368
result=dict(status='passed_all_residual_reference_summary_arithmetic',responses=15984,diagnostic_dispositions=len(frame),summary_rows=len(table),
    unresolved=int((frame.status!='descriptive_marginal_residual_diagnostics').sum()),
    source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),
    scope='Every serialized replicate count and summary mean/quantile checked; source response replay was performed by producer. No independent covariance replay or fungal-design calibration.')
write_json(Path('metadata/matched_simulation_residual_reference_readback_20260928.json'),result)
print(json.dumps(result))
