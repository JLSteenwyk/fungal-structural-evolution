"""Export all conditional working-model estimates with unchanged review dispositions."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
import psutil
from readback_full_analytic_stationarity import check_row
from screen_duplication_domain_alignment_coverage import sha

COVARIATES=['identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference']


def extract(payload, assessment):
    reasons=check_row(assessment,payload)
    row=dict(original_status=payload['status'],analytic_assessment=assessment['assessment'],
        numerical_review_required=bool(reasons),review_reasons=';'.join(reasons),
        fit_error_type=payload.get('error_type',''),intercept=None,conditional_intercept_variance=None,
        records=None,profiled_scale=None,negative_profiled_reml=None)
    for name in COVARIATES:row['coefficient_'+name]=None
    for name in ['background','family_component','species']:row['variance_ratio_'+name]=None
    if payload['status']=='fit_error_requires_review':return row
    active=payload['active_covariates'];assert len(active)==4 and all(type(v) is bool for v in active)
    names=['intercept']+[n for n,a in zip(COVARIATES,active) if a]
    beta=payload['raw_unit_beta'];cov=np.asarray(payload['raw_unit_conditional_beta_covariance'])
    assert len(beta)==len(names) and cov.shape==(len(names),len(names))
    assert all(math.isfinite(v) for v in beta) and np.isfinite(cov).all() and cov[0,0]>0
    assert payload['component_order']==['background','family_component','species']
    for name,value in zip(names,beta):row[name if name=='intercept' else 'coefficient_'+name]=value
    for name,value in zip(payload['component_order'],payload['ratios']):row['variance_ratio_'+name]=value
    row.update(conditional_intercept_variance=float(cov[0,0]),records=payload['records'],profiled_scale=payload['profiled_scale'],negative_profiled_reml=payload['negative_profiled_reml'])
    return row


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    config=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in config['pins'].items():assert sha(path)==h,path
    verify();launch=json.loads(Path(config['dependency_launch']).read_text())
    while True:
        try:
            process=psutil.Process(launch['pid'])
            if process.create_time()!=launch['created'] or process.status()==psutil.STATUS_ZOMBIE:break
            assert process.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    raw=subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True)
    assert dict(l.split('=',1) for l in raw.splitlines())==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    verify();readback=Path(config['analytic_readback']);rr=json.loads((readback/'receipt.json').read_text())
    assert rr['status']=='passed_full_analytic_stationarity_output_readback' and rr['dispositions']==144040
    source=Path(config['analytic_source']);sr=json.loads((source/'receipt.json').read_text())
    assert rr['source_receipt_sha256']==sha(source/'receipt.json')
    root=Path(config['production']);pr=json.loads((root/'receipt.json').read_text())
    assert sr['production_receipt_sha256']==sha(root/'receipt.json')
    bindings={str(folder/'receipt.json'):sha(folder/'receipt.json') for folder in [readback,source,root]}
    for folder,r in [(readback,rr),(source,sr),(root,pr)]:
        for name,h in r['artifacts'].items():assert sha(folder/name)==h;bindings[str(folder/name)]=h
    assessments={}
    for line in (source/'analytic_stationarity.jsonl').open():
        r=json.loads(line);key=r['fit_input_id'],r['tree'];assert key not in assessments;assessments[key]=r
    rows=[];seen=set()
    for line in (root/'fit_manifest.jsonl').open():
        entry=json.loads(line);key=entry['fit_input_id'],entry['tree'];assert key in assessments and key not in seen;seen.add(key)
        assert sha(entry['path'])==entry['sha256']==assessments[key]['source_fit_sha256']
        saved=json.loads(Path(entry['path']).read_text());assert (saved['fit_input_id'],saved['tree'])==key
        row=extract(saved['payload'],assessments[key]);row.update(fit_input_id=key[0],tree=key[1],source_fit_sha256=entry['sha256'])
        assert row['original_status']==entry['status'];rows.append(row)
        if len(rows)%10000==0:print('Exported unique fits',len(rows),'/ 144040',flush=True)
    assert seen==set(assessments) and len(seen)==144040
    frame=pd.DataFrame(rows)
    assert int(frame.numerical_review_required.sum())==rr['remaining_review_cases']
    support=Path(config['support']);mapping=pd.read_csv(support/'full_setting_projection.tsv',sep='\t')
    assert len(mapping)==82944 and mapping.fit_input_id.nunique()==28808
    trees=set(frame.tree);assert len(trees)==5
    assert seen=={(i,t) for i in set(mapping.fit_input_id) for t in trees}
    # Records are inherited from the inventory for failed fits; never discard failures.
    full=mapping.merge(frame,on='fit_input_id',how='left',validate='many_to_many',suffixes=('_inventory','_fit'))
    assert len(full)==414720
    has=full.records_fit.notna();assert (full.loc[has,'records_fit']==full.loc[has,'records_inventory']).all()
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    for name,table in [('unique_fits.parquet',frame),('full_settings.parquet',full)]:
        table.to_parquet(out/name,index=False);pd.testing.assert_frame_equal(pd.read_parquet(out/name),table)
    # Every expanded numeric/status field must be constant within its unique fit.
    for column in frame.columns:
        if column in ['fit_input_id','tree']:continue
        target='records_fit' if column=='records' else column
        assert full.groupby(['fit_input_id','tree'])[target].nunique(dropna=False).eq(1).all(),column
    verify()
    for path,h in bindings.items():assert sha(path)==h,path
    result=dict(status='complete_full_working_model_grid_export',plan_sha256=ph,unique_fits=len(frame),full_setting_fits=len(full),remaining_review_unique_fits=int(frame.numerical_review_required.sum()),remaining_review_setting_fits=int(full.numerical_review_required.sum()),source_bindings=bindings,artifacts={p.name:sha(p) for p in out.iterdir()},scope='All fits and original settings retained with coefficients in original units, conditional intercept variance, variance ratios, original flags and analytic review reasons. Omitted coefficients are null, not zero. Failed fits retain inventory records. No confidence intervals, p-values, biological effect or uncertainty calibration claimed.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

if __name__=='__main__':main()
