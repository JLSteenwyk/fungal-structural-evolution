"""Full interval serialization, source mapping and scalar arithmetic readback.

Does not recompute covariance contractions or establish statistical calibration.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
import psutil
from scipy.stats import t
from ancestral_chain_attempt import sha,write_json

NAMES=['intercept','identity_difference','original_coverage_difference',
       'log_aligned_length_ratio','confidence_fraction_difference']
REVIEW={'information_requires_review','adjusted_covariance_requires_review',
        'scalar_moments_require_review','scalar_interval_requires_review'}


def check_fit(fit,row):
    for key in ['fit_input_id','tree','selected_source_sha256','selection']:
        assert fit[key]==row[key],key
    status=fit['status'];intervals=fit['intervals']
    if status in ['selected_fit_requires_review','interval_computation_requires_review']:
        assert intervals=={}
        if status=='selected_fit_requires_review':assert row['selected_review_required']
        else:assert fit['error_type'] and isinstance(fit['error'],str)
        return [dict(coefficient=name,status=status,estimate=None) for name in NAMES]
    assert status=='candidate_interval_dispositions_pending_validation'
    assert not row['selected_review_required']
    assert set(intervals)==set(NAMES)
    eigen=np.asarray(fit['information_eigenvalues'])
    assert eigen.shape==(4,) and np.isfinite(eigen).all()
    assert 0<=fit['information_rank']<=4
    assert fit['information_method'] in ['nested_low_rank_trace','blocked_reference_cancellation_fallback']
    zeros=[False]+[row['selected_variance_ratio_'+name]==0 for name in ['background','family_component','species']]
    rows=[]
    for name in NAMES:
        interval=intervals[name];kind=interval['status']
        value=row['selected_intercept' if name=='intercept' else 'selected_coefficient_'+name]
        if kind=='omitted_constant_covariate':
            assert name!='intercept' and pd.isna(value) and interval==dict(status=kind,estimate=None)
        elif kind in REVIEW:
            assert not pd.isna(value) and interval['reason']
        else:
            assert kind=='candidate_interval_pending_coverage_validation'
            values=[interval[k] for k in ['estimate','adjusted_variance','denominator_df','A2','lower','upper']]
            assert np.isfinite(values).all()
            assert interval['adjusted_variance']>0 and interval['denominator_df']>0 and interval['A2']>0
            assert interval['nominal_level']==.95 and interval['numerator_df']==1 and interval['f_scaling']==1
            np.testing.assert_allclose(interval['estimate'],value,rtol=1e-7,atol=1e-8)
            np.testing.assert_allclose(interval['denominator_df'],2/interval['A2'],rtol=1e-12,atol=1e-12)
            half=t.ppf(.975,interval['denominator_df'])*np.sqrt(interval['adjusted_variance'])
            np.testing.assert_allclose([interval['lower'],interval['upper']],
                [interval['estimate']-half,interval['estimate']+half],rtol=1e-10,atol=1e-10)
            assert interval['exact_zero_components']==zeros
        rows.append(dict(coefficient=name,**interval))
    return rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args()
    config=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for name,h in config['pins'].items():assert sha(name)==h,name
    verify();launch=json.loads(Path(config['producer_launch']).read_text())
    while True:
        state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],
            '-p','MainPID','-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:
            assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0',state
            break
        assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid']
        try:
            p=psutil.Process(launch['pid']);assert p.create_time()==launch['created'] and p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:time.sleep(1);continue
        print('Waiting for complete selected interval grid',flush=True);time.sleep(30)
    verify();plan=json.loads(Path(config['producer_plan']).read_text())
    assert launch['plan_sha256']==sha(config['producer_plan'])
    for name,h in plan['pins'].items():assert sha(name)==h
    root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='complete_candidate_interval_dispositions_pending_full_readback'
    assert receipt['plan_sha256']==sha(config['producer_plan'])
    for name,h in receipt['artifacts'].items():assert sha(root/name)==h
    prep=Path(plan['inputs']);pr=json.loads((prep/'receipt.json').read_text())
    for name,h in pr['artifacts'].items():assert sha(prep/name)==h
    for name,h in pr['source_bindings'].items():assert sha(name)==h
    tasks={}
    for line in (prep/'tasks.jsonl').open():
        task=json.loads(line);assert task['fit_input_id'] not in tasks;tasks[task['fit_input_id']]=task
    frame=pd.read_parquet(plan['selected']);assert not frame.duplicated(['fit_input_id','tree']).any()
    selected={(r['fit_input_id'],r['tree']):r for r in frame.to_dict('records')}
    records=[];seen=set();fits_seen=set();fit_counts=Counter()
    for line in (root/'manifest.jsonl').open():
        entry=json.loads(line);identifier=entry['fit_input_id'];assert identifier not in seen;seen.add(identifier)
        path=Path(entry['path']);assert path.resolve().is_relative_to(root.resolve())
        assert sha(path)==entry['sha256']==path.with_suffix('.sha256').read_text().strip()
        result=json.loads(path.read_text());assert result['binding']==dict(plan_sha256=receipt['plan_sha256'],**tasks[identifier])
        assert set(result['fits'])==set(tasks[identifier]['selected_sources'])==set(pr['factors'])
        for tree,fit in result['fits'].items():
            key=(identifier,tree);assert key not in fits_seen;fits_seen.add(key)
            source=tasks[identifier]['selected_sources'][tree];assert sha(source['path'])==source['sha256']
            rows=check_fit(fit,selected[key]);fit_counts[fit['status']]+=1
            records.extend(dict(fit_input_id=identifier,tree=tree,selection=fit['selection'],**r) for r in rows)
    assert len(seen)==28808 and seen==set(tasks) and fits_seen==set(selected) and len(fits_seen)==144040
    assert len(records)==720200
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    table=pd.DataFrame(records);table.to_parquet(out/'coefficient_dispositions.parquet',index=False)
    verify();write_json(out/'receipt.json',dict(status='complete_selected_interval_mapping_and_scalar_arithmetic_readback',
        plan_sha256=ph,producer_receipt_sha256=sha(root/'receipt.json'),inputs=len(seen),fits=len(fits_seen),
        coefficient_dispositions=len(records),fit_status_counts=dict(fit_counts),coefficient_status_counts=table.status.value_counts().to_dict(),
        artifacts={'coefficient_dispositions.parquet':sha(out/'coefficient_dispositions.parquet')},
        scope='Full source/mapping/completeness and scalar arithmetic checks. Covariance contractions are not recomputed here; '
              'no coverage qualification, model adequacy or multiple-testing inference established.'))

if __name__=='__main__':main()
