"""Wait for fixed-size simulation completion, then replay all saved dispositions."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import psutil
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha,write_json
from simulate_matched_working_model import simulate
from refit_matched_simulation import response_seed
from matched_mixed_covariance import MatchedCovariance,profiled_reml
from calibrate_matched_kr import evaluate_refit,summarize
from matched_calibration_intervals import summarize_replicates


def same(actual, expected):
    if isinstance(expected,dict):
        assert actual.keys()==expected.keys()
        for key in expected:same(actual[key],expected[key])
    elif isinstance(expected,list):
        assert len(actual)==len(expected)
        for a,b in zip(actual,expected):same(a,b)
    elif isinstance(expected,float):
        np.testing.assert_allclose(actual,expected,rtol=1e-8,atol=1e-9)
    else:assert actual==expected


def replay(saved, arrays, case, plan, ph):
    assert saved['plan_sha256']==ph and saved['case']==case
    r=saved['refit'];rep=saved['replicate']
    assert r['fit_id']==case['id'] and r['replicate']==rep
    entropy=response_seed(plan['master_seed'],case['id'],rep)
    assert r['seed_entropy']==entropy
    a=arrays
    y=simulate(a['background'],a['family'],a['factor'],a['design'],a['beta'],case['scale'],case['ratios'],
               1,np.random.default_rng(np.random.SeedSequence(entropy)))[:,0]
    assert hashlib.sha256(np.asarray(y,dtype='<f8').tobytes()).hexdigest()==r['response_sha256']
    if r['status']!='refit_error_requires_review':
        f=r['fit'];assert len(f['candidates'])==24
        best=min(f['candidates'],key=lambda c:(c['objective'],not c['success']))
        np.testing.assert_array_equal(f['log1p_ratios'],best['theta'])
        value=profiled_reml(MatchedCovariance(a['background'],a['family'],a['factor'],1.,*f['ratios']),a['design'],y)
        for key in ['beta','profiled_scale','conditional_beta_covariance','negative_profiled_reml']:
            np.testing.assert_allclose(value[key],f[key],rtol=1e-7,atol=1e-8)
        if r['status']=='refit_numerically_checked':
            from scipy.stats import t
            z=(value['beta']-a['beta'])/np.sqrt(np.diag(value['conditional_beta_covariance']))
            np.testing.assert_array_equal(abs(z)<=t.ppf(.975,value['residual_degrees_of_freedom']),
                                           r['nominal_095_t_interval_covers_truth'])
    recomputed=evaluate_refit(r,a['background'],a['family'],a['factor'],a['design'],a['beta'])
    same(recomputed,saved['interval'])
    return r,recomputed


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();config=json.loads(args.plan.read_text());ch=sha(args.plan)
    threadpool_limits(1)
    def verify():
        assert sha(args.plan)==ch
        for path,digest in config['pins'].items():assert sha(path)==digest,path
    verify();launch=json.loads(Path(config['producer_launch']).read_text())
    while True:
        state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],
            '-p','ActiveState','-p','MainPID','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:break
        assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid']
        try:
            process=psutil.Process(launch['pid'])
            assert process.create_time()==launch['created'] and process.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:
            time.sleep(1);continue
        time.sleep(30)
    assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0',state
    verify();plan=json.loads(Path(config['producer_plan']).read_text());ph=sha(config['producer_plan'])
    assert ph==launch['plan_sha256']
    for path,digest in plan['pins'].items():assert sha(path)==digest,path
    root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='completed_fixed_size_synthetic_simulation_dispositions_pending_audit'
    assert receipt['plan_sha256']==ph
    for name,digest in receipt['artifacts'].items():assert sha(root/name)==digest,name
    rows=json.loads((root/'manifest.json').read_text())
    expected={(c['id'],r) for c in plan['cases'] for r in range(plan['replicates'])}
    assert len(rows)==len(expected)==receipt['attempted']==15984
    assert {(r['case'],r['replicate']) for r in rows}==expected
    summaries={};baseline={};counts=Counter();n=0
    for case in plan['cases']:
        with np.load(case['design'],allow_pickle=False) as h:a={k:h[k] for k in h.files}
        candidates=[];refits=[]
        for row in [r for r in rows if r['case']==case['id']]:
            path=Path(row['path'])
            assert path.resolve().is_relative_to(root.resolve())
            assert sha(path)==row['sha256']==path.with_suffix('.sha256').read_text().strip()
            saved=json.loads(path.read_text());assert saved['replicate']==row['replicate']
            r,interval=replay(saved,a,case,plan,ph)
            assert row['refit_status']==r['status'] and row['interval_status']==interval['status']
            refits.append(r);candidates.append(interval);counts[r['status']]+=1;n+=1
            if n%500==0:print('Audited simulated refits',n,'/15984',flush=True)
        summaries[case['id']]=summarize(candidates)
        baseline[case['id']]=summarize_replicates(refits,case['coefficients'])
    same(summaries,json.loads((root/'coverage_summaries.json').read_text()))
    verify()
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    write_json(out/'kr_coverage.json',summaries);write_json(out/'conditional_t_coverage.json',baseline)
    write_json(out/'receipt.json',dict(status='completed_all_synthetic_response_fit_interval_replays',
        plan_sha256=ch,source_receipt_sha256=sha(root/'receipt.json'),attempted=n,status_counts=dict(counts),
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='All15984 seed/response hashes, saved-source bindings, selected-fit likelihoods and '
              'candidate intervals replayed; failed outcomes retained. No independent refitting '
              'or all-candidate gradient reoptimization. Synthetic marginal coverage only; '
              'not full-real-grid calibration or model-adequacy evidence.'))


if __name__=='__main__':main()
