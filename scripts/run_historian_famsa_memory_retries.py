#!/usr/bin/env python3
"""Run full-family Historian computational sensitivities, not fitted ancestors."""
import concurrent.futures, fcntl, json, math, os, re, signal, subprocess, time
from pathlib import Path
import psutil
from Bio import Phylo, SeqIO
from prepare_case_ancestral_neighborhoods import sha


def validate(job, path):
    data = json.loads(path.read_text())
    tree = Phylo.read(job['tree'], 'newick')
    expected = {(p.name, c.name): c.branch_length for p in tree.find_clades() for c in p.clades}
    branches = {(p, c): v for p, c, v in data['branches']}
    assert len(branches) == len(data['branches']) and set(branches) == set(expected)
    assert data['root'] == job['root'] == tree.root.name
    errors = [abs(branches[k]-v) for k, v in expected.items()]
    assert all(math.isfinite(branches[k]) and abs(branches[k]-v) <= 1e-5*max(1.,abs(v)) for k,v in expected.items())
    rows = data['rowData']; assert set(rows) == {n.name for n in tree.find_clades()}
    assert all(isinstance(s,str) for s in rows.values())
    assert len({len(s) for s in rows.values()}) == 1
    records = list(SeqIO.parse(job['alignment'], 'fasta'))
    assert len(records) == job['proteins']
    imputed = []
    for r in records:
        original = str(r.seq).replace('-','').upper()
        result = rows[r.id].replace('-','').upper()
        assert len(original) == len(result), r.id
        for i,(a,b) in enumerate(zip(original,result)):
            if a == 'X': imputed.append(dict(protein=r.id,position=i+1,output=b))
            else: assert a == b, (r.id,i,a,b)
    return dict(tips_verified=len(records), output_columns=len(next(iter(rows.values()))),
                maximum_serialized_branch_error=max(errors), unknown_tip_positions=imputed)


def main():
    pp=Path('metadata/historian_famsa_memory_retry_plan_20260927.json');plan=json.loads(pp.read_text());digest=sha(pp)
    for p,h in plan['pins'].items(): assert sha(p)==h,p
    source=json.loads(Path(plan['input_receipt']).read_text())
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True)
    lock=(out/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    def run(job):
        folder=out/job['job_id'];rp=folder/'receipt.json'
        if rp.exists():
            r=json.loads(rp.read_text());assert r['plan_sha256']==digest
            for p,h in r['artifacts'].items(): assert sha(folder/p)==h
            return r
        folder.mkdir(exist_ok=False)
        assert sha(job['tree'])==job['tree_sha256'] and sha(job['alignment'])==job['alignment_sha256']
        command=[plan['binary'],'recon','-guide',str(Path(job['alignment']).resolve()),'-tree',str(Path(job['tree']).resolve())]+plan['options']+['-savemodel',str((folder/'fixed-model.json').resolve())]
        (folder/'command.json').write_text(json.dumps(command,indent=2)+'\n')
        start=time.monotonic();peak=0;status='pending_validation';details={}
        with (folder/'reconstruction.json').open('w') as stdout, (folder/'run.log').open('w') as stderr:
            proc=subprocess.Popen(command,stdout=stdout,stderr=stderr,start_new_session=True)
            identity=psutil.Process(proc.pid)
            (folder/'process.json').write_text(json.dumps(dict(pid=proc.pid,create_time=identity.create_time(),command=command),indent=2)+'\n')
            while proc.poll() is None:
                try: rss=sum(p.memory_info().rss for p in [identity]+identity.children(recursive=True) if p.is_running())
                except psutil.NoSuchProcess: rss=0
                peak=max(peak,rss)
                if time.monotonic()-start>plan['timeout_seconds'] or rss>plan['per_job_memory_bytes']:
                    status='capacity_timeout' if time.monotonic()-start>plan['timeout_seconds'] else 'capacity_memory_limit'
                    try: os.killpg(proc.pid,signal.SIGKILL)
                    except ProcessLookupError: pass
                    break
                time.sleep(1)
            code=proc.wait()
        if status=='pending_validation':
            if code!=0: status='program_failure'
            else:
                try: details=validate(job,folder/'reconstruction.json');status='input_preservation_verified_fixed_parameter_diagnostic'
                except (AssertionError,ValueError,KeyError,TypeError) as error:
                    status='output_validation_failure';details={'error':repr(error)}
        log=(folder/'run.log').read_text(errors='replace')
        details['guide_relaxation_messages']=[line for line in log.splitlines() if 'Zero forward likelihood' in line]
        details['final_likelihood_messages']=[line for line in log.splitlines() if 'Final Forward log-likelihood' in line]
        r=dict(job=job,status=status,exit_code=code,elapsed_seconds=time.monotonic()-start,peak_sampled_rss_bytes=peak,plan_sha256=digest,details=details,artifacts={p.name:sha(p) for p in folder.iterdir() if p.is_file()})
        rp.write_text(json.dumps(r,indent=2)+'\n');print(job['job_id'],status,flush=True);return r
    with concurrent.futures.ThreadPoolExecutor(max_workers=plan['workers']) as pool:results=list(pool.map(run,source['jobs']))
    counts={s:sum(r['status']==s for r in results) for s in {r['status'] for r in results}}
    receipt=dict(status='frozen_memory_retries_finished_review_required',jobs=len(results),status_counts=counts,plan_sha256=digest,job_receipts={str(p.relative_to(out)):sha(p) for p in out.glob('*/receipt.json')},scope=plan['scope'])
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
