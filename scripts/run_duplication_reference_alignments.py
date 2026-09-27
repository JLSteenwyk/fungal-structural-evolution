#!/usr/bin/env python3
"""Checkpoint additional reference pairs using two separately audited input collections."""
import argparse,csv,fcntl,hashlib,itertools,json,shutil,subprocess,time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import psutil
from run_cross_clan_alignments import parse_output
from run_ortholog_pair_guide_comparison import sha


from run_duplication_alignments import run_job
from duplication_reference_alignment_handoff import load_handoff


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Plan changed')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed pin '+path)
    verify()
    for spec in plan['input_sources'].values():
        dep=spec['producer']
        while psutil.pid_exists(dep['pid']):
            try:
                proc=psutil.Process(dep['pid'])
                if abs(proc.create_time()-dep['created'])>.01 or proc.status()==psutil.STATUS_ZOMBIE:break
                if proc.cmdline()!=dep['cmdline']:raise ValueError('Dependency identity changed')
            except psutil.NoSuchProcess:break
            time.sleep(30)
    verify();inputs,pairs,bindings,mh=load_handoff(plan)
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True);lock=(out/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (out/'receipt.json').exists():raise FileExistsError('Completed run requires review before restart')
    def jobs():
        for row in pairs:
            a=(row['model_a'],int(row['version_a']));b=(row['model_b'],int(row['version_b']))
            for mask in ['full','plddt70']:
                for order,(left,right) in enumerate([(a,b),(b,a)]):yield(row['pair_key'],left,right,mask,order)
    stream=iter(jobs());totals=Counter();completed=0
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool,(out/'checkpoint_manifest.tsv').open('w') as log:
        writer=csv.writer(log,delimiter='\t',lineterminator='\n');writer.writerow(['path','sha256','status'])
        while batch:=list(itertools.islice(stream,64)):
            if shutil.disk_usage(out).free<plan['minimum_free_disk_gib']*2**30:raise RuntimeError('Disk reserve reached')
            futures=[pool.submit(run_job,job,inputs,plan,ph,mh) for job in batch]
            for job,future in zip(batch,futures):
                path,status=future.result();totals[job[3]+':'+status]+=1;completed+=1;writer.writerow([str(path.relative_to(out)),sha(path),status])
            log.flush();(out/'state.json').write_text(json.dumps(dict(stage='aligning',completed=completed,total=4*len(pairs),counts=dict(totals)))+'\n')
            print('Directed dispositions',completed,'/',4*len(pairs),flush=True)
    verify()
    for path,digest in bindings.items():
        if sha(path)!=digest:raise ValueError('Materialization changed during alignments')
    result=dict(status='complete_reference_alignment_dispositions_pending_readback',plan_sha256=ph,input_bindings=bindings,input_bundle_sha256=mh,distinct_model_pairs=len(pairs),directed_dispositions=completed,counts=dict(totals),artifacts={'checkpoint_manifest.tsv':sha(out/'checkpoint_manifest.tsv')},scope='Both input orders of every additional reference model pair, full and pLDDT70 masks. Existing primary pairs are excluded from computation and must be joined from the primary results later; event links remain in the reference inventory. Missing/short inputs, native errors, parse errors and timeouts retained without automatic substitution/retry. Raw native output and per-job timings saved. Independent numeric readback, confidence/coverage interpretation, matched controls, domain/orientation checks and biological duplication/asymmetry tests remain pending.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
