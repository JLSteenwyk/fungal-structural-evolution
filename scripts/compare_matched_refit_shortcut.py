"""Repeat pinned refit timing responses with the exact-gradient shortcut."""
import argparse
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
from refine_matched_reml_analytic_fast import refine


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    config=json.loads(args.plan.read_text());ch=sha(args.plan);threadpool_limits(1)
    def verify():
        assert sha(args.plan)==ch
        for path,h in config['pins'].items():assert sha(path)==h,path
    verify();launch=json.loads(Path(config['producer_launch']).read_text())
    while True:
        state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','MainPID','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:break
        assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid']
        try:
            p=psutil.Process(launch['pid']);assert p.create_time()==launch['created'] and p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:time.sleep(1);continue
        time.sleep(30)
    assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0'
    plan=json.loads(Path(config['producer_plan']).read_text());source=Path(plan['output']);receipt=json.loads((source/'receipt.json').read_text())
    assert receipt['status']=='complete_pinned_design_refit_timing_dispositions' and receipt['plan_sha256']==sha(config['producer_plan'])
    assert len(receipt['jobs'])==len(plan['jobs'])==15
    for path,h in plan['pins'].items():assert sha(path)==h,path
    expected={(r['fit_input_id'],r['tree']):r for r in plan['jobs']};seen=set()
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=False);rows=[]
    for record in receipt['jobs']:
        key=record['fit_input_id'],record['tree'];assert key not in seen;seen.add(key)
        path=source/record['path'];assert sha(path)==record['sha256'];prior=json.loads(path.read_text());job=prior['job'];assert job==expected[key]
        assert prior['result']['status']=='refit_numerically_checked'
        with np.load(job['cache'],allow_pickle=False) as h:a={k:h[k] for k in h.files}
        with np.load(job['factor'],allow_pickle=False) as h:factor=h['factor'][a['pattern_rows']]
        original=json.loads(Path(job['original_fit']).read_text())['payload']
        x=np.column_stack([np.ones(len(a['matrix'])),a['matrix'][:,1:][:,a['active_covariates']]/a['covariate_scales']])
        started=time.perf_counter()
        response=simulate(a['background'],a['family'],factor,x,original['beta'],original['profiled_scale'],original['ratios'],1,
                          np.random.default_rng(np.random.SeedSequence(prior['result']['seed_entropy'])))[:,0]
        assert hashlib.sha256(np.asarray(response,dtype='<f8').tobytes()).hexdigest()==prior['result']['response_sha256']
        fitted=refine(a['background'],a['family'],factor,x,response,np.log1p(original['ratios']))
        elapsed=time.perf_counter()-started;old=prior['result']['fit']
        assert [(c['active'],c['start']) for c in fitted['candidates']]==[(c['active'],c['start']) for c in old['candidates']]
        checks=dict(objective_agrees=bool(abs(fitted['negative_profiled_reml']-old['negative_profiled_reml'])<=1e-6),
            beta_agrees=bool(np.allclose(fitted['beta'],old['beta'],rtol=1e-5,atol=1e-7)),
            covariance_agrees=bool(np.allclose(fitted['conditional_beta_covariance'],old['conditional_beta_covariance'],rtol=1e-5,atol=1e-7)),
            status_agrees=fitted['status']==old['status'],diagnostic_flags_agree=fitted['checks']==old['checks'])
        target=out/path.name
        write_json(target,dict(plan_sha256=ch,source_sha256=sha(path),fit=fitted,checks=checks,elapsed_seconds=elapsed))
        rows.append(dict(fit_input_id=key[0],tree=key[1],records=len(x),checks=checks,
            original_seconds=prior['simulation_and_refit_seconds'],shortcut_seconds=elapsed,
            speed_ratio=prior['simulation_and_refit_seconds']/elapsed,
            objective_difference=fitted['negative_profiled_reml']-old['negative_profiled_reml'],
            path=target.name,sha256=sha(target)))
        print(len(rows),'/15',checks,'speed',round(rows[-1]['speed_ratio'],3),flush=True)
    assert seen==set(expected);verify()
    write_json(out/'receipt.json',dict(status='complete_paired_shortcut_refit_comparison',plan_sha256=ch,
        source_receipt_sha256=sha(source/'receipt.json'),jobs=rows,
        all_comparisons_pass=all(all(r['checks'].values()) for r in rows),
        scope='Same15 response hashes and24-candidate grids. Review every discrepancy; no global equivalence, representative full-grid speedup, calibration or production replacement claim.'))


if __name__=='__main__':main()
