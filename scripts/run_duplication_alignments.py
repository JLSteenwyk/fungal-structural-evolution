#!/usr/bin/env python3
"""Checkpoint both input orders and both confidence masks for every queued model pair."""
import argparse,csv,fcntl,hashlib,itertools,json,shutil,subprocess,time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import psutil
from run_cross_clan_alignments import parse_output
from run_ortholog_pair_guide_comparison import sha


def run_job(job,inputs,plan,ph,manifest_hash):
    pair,left,right,mask,order=job;rows=[inputs[(*key,mask)] for key in [left,right]]
    out=Path(plan['output'])/'pairs'/pair[:2];out.mkdir(parents=True,exist_ok=True);dest=out/f'{pair}-{mask}-{order}.json'
    identities=[dict(model_id=key[0],version=key[1],status=row['status'],path=row.get('path'),sha256=row.get('sha256')) for key,row in zip([left,right],rows)]
    ready=all(r['status']=='ready' for r in rows)
    for row in rows:
        if row['status']=='ready' and sha(row['path'])!=row['sha256']:raise ValueError('Changed alignment coordinate')
    command=[plan['usalign'],*[r['path'] for r in rows],*plan['options']] if ready else None
    binding=dict(plan_sha256=ph,input_manifest_sha256=manifest_hash,pair_key=pair,mask=mask,order=order,inputs=identities,command=command)
    if dest.exists():
        record=json.loads(dest.read_text())
        if any(record[k]!=v for k,v in binding.items()):raise ValueError('Checkpoint binding differs')
        if record['status']=='aligned':
            if not ready or record['returncode']!=0 or parse_output(record['stdout'],[r['sequence'] for r in rows])!=record['metrics']:raise ValueError('Checkpoint metrics differ')
        elif record['status']=='input_unavailable':
            if ready:raise ValueError('Unexpected exclusion checkpoint')
        elif record['status'] not in ['native_error','timeout','parse_error']:raise ValueError('Unknown checkpoint disposition')
        return dest,record['status']
    record=dict(binding);started=time.monotonic()
    if not ready:record.update(status='input_unavailable',input_statuses=[r['status'] for r in rows])
    else:
        try:
            native=subprocess.run(command,capture_output=True,text=True,timeout=plan['per_pair_timeout_seconds'])
            record.update(stdout=native.stdout,stderr=native.stderr,returncode=native.returncode)
            if native.returncode:record['status']='native_error'
            else:
                try:record.update(status='aligned',metrics=parse_output(native.stdout,[r['sequence'] for r in rows]))
                except (ValueError,AttributeError,StopIteration,IndexError) as e:record.update(status='parse_error',error=f'{type(e).__name__}: {e}')
        except subprocess.TimeoutExpired as e:
            def text(x):return x.decode(errors='replace') if isinstance(x,bytes) else (x or '')
            record.update(status='timeout',stdout=text(e.stdout),stderr=text(e.stderr),timeout_seconds=plan['per_pair_timeout_seconds'])
    record['elapsed_seconds']=time.monotonic()-started
    for row in rows:
        if row['status']=='ready' and sha(row['path'])!=row['sha256']:raise ValueError('Coordinate changed during alignment')
    temp=dest.with_suffix('.partial')
    if temp.exists():raise FileExistsError('Unreceipted checkpoint requires review')
    temp.write_text(json.dumps(record,separators=(',',':'))+'\n');temp.replace(dest)
    return dest,record['status']


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Plan changed')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed pin '+path)
    verify();dep=plan['producer']
    while psutil.pid_exists(dep['pid']):
        try:
            proc=psutil.Process(dep['pid'])
            if abs(proc.create_time()-dep['created'])>.01 or proc.status()==psutil.STATUS_ZOMBIE:break
            if proc.cmdline()!=dep['cmdline']:raise ValueError('Dependency identity changed')
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify();folder=Path(plan['inputs']);rp=folder/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    if r['status']!='complete_duplication_alignment_input_materialization' or r['plan_sha256']!=sha(plan['input_plan']):raise ValueError('Incomplete materialization')
    queue=Path(plan['queue']);qr=json.loads((queue/'receipt.json').read_text())
    if r['queue_receipt_sha256']!=sha(queue/'receipt.json'):raise ValueError('Queue differs from materialization')
    pairsfile=queue/'model_pairs.tsv';modelsfile=queue/'models.jsonl'
    for path in [pairsfile,modelsfile]:
        if sha(path)!=qr['artifacts'][path.name]:raise ValueError('Changed queue artifact')
    with modelsfile.open() as f:expected={(m['model_id'],m['version'],mask) for line in f for m in [json.loads(line)] for mask in ['full','plddt70']}
    manifest=folder/'inputs.jsonl';mh=sha(manifest)
    if mh!=r['artifacts']['inputs.jsonl']:raise ValueError('Changed input manifest')
    inputs={};counts=Counter()
    with manifest.open() as f:
        for line in f:
            row=json.loads(line);key=(row['model_id'],row['version'],row['mask'])
            if key in inputs or row['status'] not in ['ready','too_few_retained_residues','source_rejected']:raise ValueError('Invalid input disposition')
            counts[row['mask']+':'+row['status']]+=1
            inputs[key]={k:row[k] for k in ['status','path','sha256','sequence'] if k in row}
    if set(inputs)!=expected or dict(counts)!=r['counts']:raise ValueError('Input universe differs')
    with pairsfile.open() as f:pairs=list(csv.DictReader(f,delimiter='\t'))
    seen=set()
    for row in pairs:
        endpoints=[(row['model_a'],int(row['version_a'])),(row['model_b'],int(row['version_b']))]
        key=hashlib.sha256(json.dumps(sorted(endpoints),separators=(',',':')).encode()).hexdigest()
        if key!=row['pair_key'] or key in seen or endpoints[0]==endpoints[1]:raise ValueError('Invalid pair identity')
        for endpoint in endpoints:
            for mask in ['full','plddt70']:
                if (*endpoint,mask) not in inputs:raise ValueError('Missing pair input')
        seen.add(key)
    if len(pairs)!=qr['unique_distinct_model_pairs']:raise ValueError('Pair count differs')
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
    if sha(manifest)!=mh or sha(rp)!=rh:raise ValueError('Materialization changed during alignments')
    result=dict(status='complete_duplication_alignment_dispositions_pending_readback',plan_sha256=ph,input_receipt_sha256=rh,input_manifest_sha256=mh,distinct_model_pairs=len(pairs),directed_dispositions=completed,counts=dict(totals),artifacts={'checkpoint_manifest.tsv':sha(out/'checkpoint_manifest.tsv')},scope='Both input orders of every distinct model pair, full and pLDDT70 masks. Missing/short inputs, native errors, parse errors and timeouts retained without automatic substitution/retry. Raw native output and per-job timings saved. Independent numeric readback, confidence/coverage interpretation, matched controls, domain/orientation checks and biological duplication/asymmetry tests remain pending.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
