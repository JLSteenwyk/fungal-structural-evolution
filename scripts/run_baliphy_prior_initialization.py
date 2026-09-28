#!/usr/bin/env python3
"""Initialize explicit priors across all effective inputs; no posterior MCMC."""
import concurrent.futures, fcntl, json, os, signal, subprocess, time
from pathlib import Path
import psutil
from prepare_case_ancestral_neighborhoods import sha


def main():
    pp=Path('metadata/baliphy_prior_initialization_plan_20260927.json');plan=json.loads(pp.read_text());digest=sha(pp)
    for p,h in plan['pins'].items():assert sha(p)==h,p
    jobs=plan['jobs'];assert len(jobs)==405
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True)
    lock=(out/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    def run(job):
        folder=out/job['job_id'];rp=folder/'receipt.json'
        if rp.exists():
            r=json.loads(rp.read_text());assert r['plan_sha256']==digest
            for f,h in r['artifacts'].items():assert sha(folder/f)==h
            return r
        folder.mkdir(exist_ok=False)
        assert sha(job['alignment'])==job['alignment_sha256'] and sha(job['tree'])==job['tree_sha256']
        command=[plan['binary'],str(Path(job['alignment']).resolve()),'--fix','tree='+str(Path(job['tree']).resolve()),'--name','initialization']+plan['options']+['--smodel',job['smodel'],'--imodel',job['imodel']]
        (folder/'command.json').write_text(json.dumps(command,indent=2)+'\n')
        start=time.monotonic();peak=0;status=None
        with (folder/'stdout.log').open('w') as stdout,(folder/'stderr.log').open('w') as stderr:
            process=subprocess.Popen(command,cwd=folder,stdout=stdout,stderr=stderr,start_new_session=True)
            identity=psutil.Process(process.pid)
            (folder/'process.json').write_text(json.dumps(dict(pid=process.pid,created=identity.create_time(),command=command),indent=2)+'\n')
            while process.poll() is None:
                try:rss=sum(p.memory_info().rss for p in [identity]+identity.children(recursive=True) if p.is_running())
                except psutil.NoSuchProcess:rss=0
                peak=max(peak,rss)
                if time.monotonic()-start>plan['timeout_seconds'] or rss>plan['memory_limit_bytes']:
                    status='initialization_timeout' if time.monotonic()-start>plan['timeout_seconds'] else 'initialization_memory_limit'
                    try:os.killpg(process.pid,signal.SIGKILL)
                    except ProcessLookupError:pass
                    break
                time.sleep(1)
            code=process.wait()
        if status is None:status='initialization_exited_zero_pending_readback' if code==0 else 'initialization_failed'
        r=dict(status=status,job=job,plan_sha256=digest,exit_code=code,elapsed_seconds=time.monotonic()-start,peak_sampled_rss_bytes=peak,artifacts={str(p.relative_to(folder)):sha(p) for p in folder.rglob('*') if p.is_file()})
        rp.write_text(json.dumps(r,indent=2)+'\n');print(job['job_id'],status,flush=True);return r
    with concurrent.futures.ThreadPoolExecutor(max_workers=plan['resources']['workers']) as pool:results=list(pool.map(run,jobs))
    counts={s:sum(r['status']==s for r in results) for s in {r['status'] for r in results}}
    receipt=dict(status='full_initialization_grid_finished_pending_readback',jobs=len(results),status_counts=counts,plan_sha256=digest,job_receipts={str(p.relative_to(out)):sha(p) for p in out.glob('*/receipt.json')},scope=plan['scope'])
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
