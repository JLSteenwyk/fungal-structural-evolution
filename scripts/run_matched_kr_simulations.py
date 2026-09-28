"""Fixed-size synthetic calibration, with all refit/interval failures retained."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import fcntl
import json
from pathlib import Path
import shutil
import time
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha, write_json
from refit_matched_simulation import refit_one
from calibrate_matched_kr import evaluate_refit, summarize


def run_job(job):
    threadpool_limits(1)
    case, replicate, plan, ph = job
    path = Path(plan['output'])/'replicates'/case['id']/f'{replicate:05d}.json'
    digest_path = path.with_suffix('.sha256')
    if path.exists() or digest_path.exists():
        assert path.exists() and digest_path.exists(), 'Incomplete saved disposition requires review'
        assert sha(path) == digest_path.read_text().strip()
        saved = json.loads(path.read_text())
        assert saved['plan_sha256'] == ph and saved['case'] == case and saved['replicate'] == replicate
    else:
        with np.load(case['design'], allow_pickle=False) as h:
            a = {k:h[k] for k in h.files}
        started = time.perf_counter()
        refit = refit_one(a['background'],a['family'],a['factor'],a['design'],a['beta'],
                          case['scale'],case['ratios'],plan['master_seed'],case['id'],replicate)
        interval = evaluate_refit(refit,a['background'],a['family'],a['factor'],a['design'],a['beta'])
        saved = dict(plan_sha256=ph,case=case,replicate=replicate,refit=refit,interval=interval,
                     elapsed_seconds=time.perf_counter()-started)
        write_json(path,saved)
        digest_path.write_text(sha(path)+'\n')
    return dict(case=case['id'],replicate=replicate,path=str(path),sha256=sha(path),
                bytes=path.stat().st_size,refit_status=saved['refit']['status'],
                interval_status=saved['interval']['status'],elapsed_seconds=saved['elapsed_seconds'])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    threadpool_limits(1)
    def verify():
        assert sha(args.plan)==ph
        for name,digest in plan['pins'].items():assert sha(name)==digest,name
    verify()
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        binding=out/'run_plan.json'
        if binding.exists():assert sha(binding)==ph
        else:shutil.copyfile(args.plan,binding)
        if (out/'receipt.json').exists():raise RuntimeError('Completed output already exists; audit it')
        assert shutil.disk_usage(out).free > plan['resources']['minimum_free_gib']*2**30
        for case in plan['cases']:(out/'replicates'/case['id']).mkdir(parents=True,exist_ok=True)
        jobs=[(case,r,plan,ph) for case in plan['cases'] for r in range(plan['replicates'])]
        assert len({(c['id'],r) for c,r,_,_ in jobs})==len(jobs)
        rows=[];total_bytes=0
        with ProcessPoolExecutor(max_workers=plan['resources']['workers']) as pool:
            for row in pool.map(run_job,jobs,chunksize=1):
                rows.append(row);total_bytes+=row['bytes']
                if len(rows)%100==0:
                    verify()
                    assert total_bytes < plan['resources']['storage_gib']*2**30
                    assert shutil.disk_usage(out).free > plan['resources']['minimum_free_gib']*2**30
                    print('Simulation dispositions',len(rows),'/',len(jobs),flush=True)
        summaries={}
        for case in plan['cases']:
            selected=[r for r in rows if r['case']==case['id']]
            assert len(selected)==plan['replicates']
            records=[]
            for row in selected:
                assert sha(row['path'])==row['sha256']
                records.append(json.loads(Path(row['path']).read_text())['interval'])
            summaries[case['id']]=summarize(records)
        verify()
        write_json(out/'manifest.json',rows);write_json(out/'coverage_summaries.json',summaries)
        write_json(out/'receipt.json',dict(status='completed_fixed_size_synthetic_simulation_dispositions_pending_audit',
            plan_sha256=ph,attempted=len(rows),cases=len(plan['cases']),replicates_per_case=plan['replicates'],
            artifacts={name:sha(out/name) for name in ['run_plan.json','manifest.json','coverage_summaries.json']},
            scope='Synthetic Gaussian working-model calibration; all dispositions retained. '
                  'No claim of coverage for all real designs, misspecification robustness, '
                  'multiple-testing control or completed biological analysis.'))


if __name__=='__main__':main()
