"""Time complete simulation/refit calls on pinned existing designs; not calibration."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha,write_json
from audit_matched_simulation_cache import replay
from refit_matched_simulation import refit_one


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan);threadpool_limits(1)
    def verify():
        assert sha(args.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    rows=[]
    for job in plan['jobs']:
        started=time.perf_counter()
        with np.load(job['cache'],allow_pickle=False) as h:a={k:h[k] for k in h.files}
        with np.load(job['factor'],allow_pickle=False) as h:factor=h['factor']
        old=json.loads(Path(job['original_fit']).read_text())['payload']
        readback=replay(a,old,factor)
        replay_seconds=time.perf_counter()-started
        x=np.column_stack([np.ones(len(a['matrix'])),a['matrix'][:,1:][:,a['active_covariates']]/a['covariate_scales']])
        factor=factor[a['pattern_rows']]
        started=time.perf_counter()
        result=refit_one(a['background'],a['family'],factor,x,np.asarray(old['beta']),old['profiled_scale'],old['ratios'],
                         20260928,job['fit_input_id']+'/'+job['tree']+'/benchmark',0)
        elapsed=time.perf_counter()-started
        path=out/(job['fit_input_id']+'-'+job['tree']+'.json')
        write_json(path,dict(job=job,plan_sha256=ph,original_replay=readback,original_replay_and_load_seconds=replay_seconds,
                            simulation_and_refit_seconds=elapsed,result=result))
        rows.append(dict(fit_input_id=job['fit_input_id'],tree=job['tree'],records=len(x),status=result['status'],
                         simulation_and_refit_seconds=elapsed,original_fit_seconds=old['wall_seconds'],path=path.name,sha256=sha(path)))
        print(len(rows),'/',len(plan['jobs']),result['status'],round(elapsed,3),'seconds',flush=True)
    verify();write_json(out/'receipt.json',dict(status='complete_pinned_design_refit_timing_dispositions',plan_sha256=ph,
        jobs=rows,scope='One simulated response for each size/tree timing case; not a coverage estimate, representative full-grid speedup or replacement for full analysis. Generation uses original fitted parameters; refined estimates are separately pending.'))


if __name__=='__main__':main()
