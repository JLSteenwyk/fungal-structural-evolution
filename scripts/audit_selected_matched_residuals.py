"""Replay all serialized residual diagnostics, retaining every disposition."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha,write_json
from run_selected_matched_residuals import wait_for_replay
from evaluate_selected_matched_residuals import evaluate


def compare(actual,expected):
    if isinstance(expected,dict):
        assert actual.keys()==expected.keys()
        for key in expected:compare(actual[key],expected[key])
    elif isinstance(expected,list):
        assert len(actual)==len(expected)
        for a,b in zip(actual,expected):compare(a,b)
    elif isinstance(expected,float):
        assert np.isfinite(expected) and np.isfinite(actual)
        np.testing.assert_allclose(actual,expected,rtol=1e-8,atol=1e-9)
    else:assert actual==expected


def replay(row,arrays,factor):
    try:return evaluate(row,arrays,factor)
    except (AssertionError,ValueError,FloatingPointError,np.linalg.LinAlgError) as error:
        return dict(**{k:row[k] for k in ['fit_input_id','tree','selected_source_sha256','selection']},
            status='residual_computation_requires_review',error_type=type(error).__name__,error=str(error))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',required=True,type=Path);a=ap.parse_args()
    config=json.loads(a.plan.read_text());ph=sha(a.plan);threadpool_limits(1)
    def verify():
        assert sha(a.plan)==ph
        for name,h in config['pins'].items():assert sha(name)==h,name
    verify();wait_for_replay({'launch':config['producer_launch']});verify()
    plan=json.loads(Path(config['producer_plan']).read_text());producer_hash=sha(config['producer_plan'])
    launch=json.loads(Path(config['producer_launch']).read_text());assert launch['plan_sha256']==producer_hash
    for name,h in plan['pins'].items():assert sha(name)==h,name
    root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='complete_descriptive_residual_dispositions_pending_full_readback'
    assert receipt['plan_sha256']==producer_hash and receipt['inputs']==28808 and receipt['fits']==144040
    for name,h in receipt['artifacts'].items():assert sha(root/name)==h
    prep=Path(plan['inputs']);r=json.loads((prep/'receipt.json').read_text())
    for name,h in r['artifacts'].items():assert sha(prep/name)==h
    for name,h in r['source_bindings'].items():assert sha(name)==h
    tasks={}
    for line in (prep/'tasks.jsonl').open():
        task=json.loads(line);assert task['fit_input_id'] not in tasks;tasks[task['fit_input_id']]=task
    factors={}
    for tree,spec in r['factors'].items():
        assert sha(spec['path'])==spec['sha256']
        with np.load(spec['path'],allow_pickle=False) as h:factors[tree]=h['factor']
    frame=pd.read_parquet(plan['selected']);assert not frame.duplicated(['fit_input_id','tree']).any()
    rows={(r['fit_input_id'],r['tree']):r for r in frame.to_dict('records')}
    seen=set();fits_seen=set();summaries=[]
    for line in (root/'manifest.jsonl').open():
        item=json.loads(line);identifier=item['fit_input_id'];assert identifier not in seen;seen.add(identifier)
        target=Path(item['path']);assert target.resolve().is_relative_to(root.resolve())
        assert sha(target)==item['sha256']==target.with_suffix('.sha256').read_text().strip()
        saved=json.loads(target.read_text());task=tasks[identifier]
        assert saved['binding']==dict(plan_sha256=producer_hash,**task)
        assert set(saved['fits'])==set(task['selected_sources'])==set(factors)
        entry=task['cache_entry'];cache=Path(plan['cache'])/entry['path']
        assert cache.resolve().is_relative_to(Path(plan['cache']).resolve()) and sha(cache)==entry['sha256']
        with np.load(cache,allow_pickle=False) as h:arrays={k:h[k] for k in h.files}
        for key,field in [('matrix','values_sha256'),('row_identity','ordered_identity_sha256')]:
            assert hashlib.sha256(arrays[key].tobytes()).hexdigest()==entry['recipe'][field]
        for tree,factor in factors.items():
            key=(identifier,tree);assert key not in fits_seen;fits_seen.add(key)
            source=task['selected_sources'][tree];assert sha(source['path'])==source['sha256']==rows[key]['selected_source_sha256']
            expected=replay(rows[key],arrays,factor);compare(saved['fits'][tree],expected)
            summary={k:expected[k] for k in ['fit_input_id','tree','selection','status']}
            if expected['status']=='descriptive_marginal_residual_diagnostics':
                summary.update({k:v for k,v in expected['summaries'].items() if not isinstance(v,list)})
                summary['records']=expected['records']
            else:summary['reason']=expected.get('reason',expected.get('error','selected_fit_requires_review'))
            summaries.append(summary)
        if len(seen)%100==0:print('Replayed residual inputs',len(seen),'/28808',flush=True)
    assert len(seen)==28808 and seen==set(tasks) and fits_seen==set(rows) and len(summaries)==144040
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    table=pd.DataFrame(summaries);table.to_parquet(out/'residual_fit_summaries.parquet',index=False)
    verify();write_json(out/'receipt.json',dict(status='complete_all_selected_residual_diagnostic_replays',
        plan_sha256=ph,producer_receipt_sha256=sha(root/'receipt.json'),inputs=len(seen),fits=len(summaries),
        status_counts=table.status.value_counts().to_dict(),artifacts={'residual_fit_summaries.parquet':sha(out/'residual_fit_summaries.parquet')},
        scope='All serialized residual summaries, quantiles and bins recomputed with the same numerical evaluator; all dispositions retained. '
              'Dense fixtures provide a separate numerical oracle. This is not independent statistical model validation or a calibrated adequacy test.'))

if __name__=='__main__':main()
