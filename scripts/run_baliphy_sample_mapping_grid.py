#!/usr/bin/env python3
"""Run short full-grid chains with in-process tree export; no convergence claim."""
import concurrent.futures, fcntl, json, os, signal, subprocess, time
from pathlib import Path
import psutil
from prepare_case_ancestral_neighborhoods import sha


def main():
    pp=Path('metadata/baliphy_sample_mapping_plan_20260927.json');plan=json.loads(pp.read_text());digest=sha(pp)
    for p,h in plan['pins'].items():assert sha(p)==h,p
    inputs=json.loads(Path(plan['input_receipt']).read_text());jobs=inputs['jobs'];assert len(jobs)==324
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
        source_root=Path(plan['initialization_output'])/job['job_id']
        source_receipt=source_root/'receipt.json'
        identity=json.loads(Path(plan['initialization_launch']).read_text())
        while not source_receipt.exists():
            state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',identity['unit'],'-p','ActiveState','-p','Result','-p','MainPID'],text=True).splitlines())
            assert state['ActiveState']=='active',state
            parent=psutil.Process(identity['pid']);assert parent.create_time()==identity['created'] and parent.cmdline()==identity['cmdline']
            time.sleep(30)
        source_result=json.loads(source_receipt.read_text())
        assert source_result['job']==job
        for f,h in source_result['artifacts'].items():assert sha(source_root/f)==h
        if source_result['exit_code']!=0:
            r=dict(status='source_initialization_unsuccessful',job=job,plan_sha256=digest,source_receipt_sha256=sha(source_receipt),artifacts={})
            rp.write_text(json.dumps(r,indent=2)+'\n');return r
        code=(source_root/'BAliPhy.Main.hs').read_text()
        marker=';mcmcState <- makeMCMCState'
        assert code.count(marker)==1
        treefile=str((folder/'runtime-tree.nwk').resolve())
        insertion=';T.writeFile '+json.dumps(treefile)+' (writeNewick_rooted (addInternalLabels tree))\n'
        code=code.replace(marker,insertion+marker)
        model=folder/'MappedModel.hs';model.write_text(code)
        (folder/'source.json').write_text(json.dumps(dict(receipt=str(source_receipt),receipt_sha256=sha(source_receipt),original_program_sha256=sha(source_root/'BAliPhy.Main.hs'),mapped_program_sha256=sha(model),change='Add same-process labeled tree export immediately before MCMC state creation; no model edits.'),indent=2)+'\n')
        command=[plan['binary'],'--seed','20260928','run',str(model.resolve()),'--iterations','20','--name','mapping-check']
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
                    status='short_chain_timeout' if time.monotonic()-start>plan['timeout_seconds'] else 'short_chain_memory_limit'
                    try:os.killpg(process.pid,signal.SIGKILL)
                    except ProcessLookupError:pass
                    break
                time.sleep(1)
            code=process.wait()
        if status is None:status='short_chain_exited_zero_pending_sample_readback' if code==0 else 'short_chain_failed'
        r=dict(status=status,job=job,plan_sha256=digest,exit_code=code,elapsed_seconds=time.monotonic()-start,peak_sampled_rss_bytes=peak,artifacts={str(p.relative_to(folder)):sha(p) for p in folder.rglob('*') if p.is_file()})
        rp.write_text(json.dumps(r,indent=2)+'\n');print(job['job_id'],status,flush=True);return r
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,jobs))
    counts={s:sum(r['status']==s for r in results) for s in {r['status'] for r in results}}
    receipt=dict(status='full_short_chain_grid_finished_pending_sample_readback',jobs=len(results),status_counts=counts,plan_sha256=digest,job_receipts={str(p.relative_to(out)):sha(p) for p in out.glob('*/receipt.json')},scope=plan['scope'])
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
