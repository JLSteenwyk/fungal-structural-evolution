#!/usr/bin/env python3
"""Run complete clade discovery with isolated outputs and identity-checked receipts."""
import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time
import psutil
from audit_orthology_family_universe import groups


def sha(p):
    digest=hashlib.sha256()
    with Path(p).open('rb') as h:
        for block in iter(lambda:h.read(8*1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def save(path,obj):
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(path)


def fasta(path):
    name=None; sequence=[]
    with path.open() as h:
        for line in h:
            if line.startswith('>'):
                if name is not None:yield name,''.join(sequence)
                name=line[1:].strip();sequence=[]
            else:sequence.append(line.strip())
    if name is not None:yield name,''.join(sequence)


def bytes_used(root):
    return sum(p.stat().st_size for p in root.rglob('*') if p.is_file())


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    a=ap.parse_args();plan=json.loads(a.plan.read_text());root=Path(plan['output']).resolve()
    root.mkdir(parents=True,exist_ok=True)
    lock=(root/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    plan_hash=sha(a.plan);stamp=root/'plan_sha256.txt'
    if stamp.exists() and stamp.read_text().strip()!=plan_hash:raise ValueError('Existing output has another plan')
    stamp.write_text(plan_hash+'\n')
    for name,h in plan['pins'].items():
        if sha(name)!=h:raise ValueError('Pin mismatch: '+name)
    limits=plan['resources'];started=time.time();completed=[]
    state=dict(status='running',pid=os.getpid(),created=psutil.Process().create_time(),plan_sha256=plan_hash,completed_clades=[])
    save(root/'state.json',state)
    env=dict(os.environ,PATH=plan['environment_path'],OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')

    def guard(starting=False):
        if shutil.disk_usage(root).free<limits['minimum_free_disk_gib']*2**30:raise RuntimeError('Disk free-space gate')
        available=psutil.virtual_memory().available/2**30
        if available<(limits['minimum_available_memory_gib'] if starting else limits['stop_available_memory_gib']):raise RuntimeError('Available memory gate')
        if bytes_used(root)>limits['output_allowance_gib']*2**30:raise RuntimeError('Output allowance exceeded')

    def command(argv,folder,label):
        guard();before=time.time()
        state.update(active_step=label,active_command=argv);save(root/'state.json',state)
        with (folder/(label+'.log')).open('w') as log:
            p=subprocess.Popen(argv,cwd=plan['working_directory'],env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            state.update(child_pid=p.pid,child_created=psutil.Process(p.pid).create_time());save(root/'state.json',state)
            try:
                while True:
                    try:rc=p.wait(timeout=30);break
                    except subprocess.TimeoutExpired:
                        guard()
                        if time.time()-before>limits['step_timeout_hours']*3600:raise TimeoutError(label)
            except BaseException:
                if p.poll() is None:
                    os.killpg(p.pid,signal.SIGTERM)
                    try:p.wait(timeout=10)
                    except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
                raise
        if rc!=0:raise RuntimeError(f'{label} returned {rc}')
        state.pop('child_pid',None);state.pop('child_created',None)
        return dict(command=argv,returncode=rc,seconds=time.time()-before,log_sha256=sha(folder/(label+'.log')))

    try:
        for case in plan['cases']:
            key=case['clade_id'];folder=root/key;receipt_path=folder/'receipt.json'
            if receipt_path.exists():
                r=json.loads(receipt_path.read_text())
                if r['plan_sha256']!=plan_hash or r['status']!='complete_identity_checked_clade_discovery':raise ValueError('Invalid completed clade receipt')
                for name,h in r['artifacts'].items():
                    if sha(folder/name)!=h:raise ValueError('Changed completed clade')
                completed.append(key);continue
            if folder.exists():raise RuntimeError('Incomplete clade directory requires review: '+key)
            guard(starting=True);folder.mkdir();state.update(active_clade=key,completed_clades=completed);save(root/'state.json',state)
            inputs=Path(case['input_directory']).resolve();expected=set();records=[]
            actual={p.name for p in inputs.iterdir() if p.is_file()}
            if actual!={f['filename'] for f in case['input_files']}:raise ValueError('Input file universe differs')
            for item in case['input_files']:
                path=inputs/item['filename']
                if sha(path)!=item['sha256']:raise ValueError('Input FASTA changed')
                count=0
                for name,seq in fasta(path):
                    if name in expected or not seq or int(name.split('_')[0])!=item['species_id']:raise ValueError('Duplicate/invalid original protein identity')
                    expected.add(name);count+=1
                    if case['taxa']==1:records.append((name,seq))
                if count!=item['proteins']:raise ValueError('Input count differs')
            if len(expected)!=case['proteins'] or len(case['input_files'])!=case['taxa']:raise ValueError('Clade dimensions differ')
            base=[plan['orthofinder'],'-t',str(limits['search_threads']),'-a',str(limits['analysis_threads']),'-S','diamond','--scores-v2','--only-groups','--no-fix-files','-I','1.2']
            steps=[];self_search=None
            if case['taxa']==1:
                wd=folder/'SearchInput';wd.mkdir()
                with (wd/'Species0.fa').open('w') as h, (wd/'SequenceIDs.txt').open('w') as ids:
                    for i,(name,seq) in enumerate(records):
                        h.write(f'>0_{i}\n{seq}\n');ids.write(f'0_{i}: {name}\n')
                (wd/'SpeciesIDs.txt').write_text('0: '+case['input_files'][0]['filename']+'\n')
                db=wd/'diamondDBSpecies0';table=wd/'Blast0_0.txt'
                steps.append(command([plan['diamond'],'makedb','--ignore-warnings','--in',str(wd/'Species0.fa'),'-d',str(db)],folder,'makedb'))
                steps.append(command([plan['diamond'],'blastp','--ignore-warnings','-d',str(db),'-q',str(wd/'Species0.fa'),'-o',str(table),'--matrix','BLOSUM62','--gapopen','11','--gapextend','1','--more-sensitive','-p',str(limits['search_threads']),'--quiet','-e','0.001','--compress','1'],folder,'self-search'))
                legal={f'0_{i}' for i in range(len(records))};hits=0
                with gzip.open(str(table)+'.gz','rt') as h:
                    for line in h:
                        fields=line.split()
                        if len(fields)!=12 or fields[0] not in legal or fields[1] not in legal or any(not math.isfinite(float(x)) for x in fields[2:]):raise ValueError('Invalid self-search row')
                        hits+=1
                self_search=dict(rows=hits,sha256=sha(str(table)+'.gz'))
                steps.append(command(base+['-b',str(wd),'-n',key],folder,'orthofinder'))
                results=list((wd/'OrthoFinder').glob('Results_*'));ids_file=wd/'SequenceIDs.txt'
            else:
                steps.append(command(base+['-f',str(inputs),'-o',str(folder/'Native')],folder,'orthofinder'))
                results=list((folder/'Native').glob('Results_*'));ids_file=None
            if len(results)!=1:raise ValueError('Expected exactly one native result directory')
            native=results[0]/'WorkingDirectory';cluster_files=list(native.glob('clusters*_id_pairs.txt'))
            if len(cluster_files)!=1:raise ValueError('Missing/ambiguous native cluster output')
            if ids_file is None:ids_file=native/'SequenceIDs.txt'
            ids={}
            with ids_file.open() as h:
                for line in h:
                    local,original=line.rstrip('\n').split(': ',1)
                    if local in ids:raise ValueError('Duplicate local mapping')
                    ids[local]=original
            if set(ids.values())!=expected or len(ids)!=len(expected):raise ValueError('Native original-ID mapping differs')
            seen=set();family_count=0;singletons=0
            export=folder/'family_membership.tsv'
            with export.open('w') as h:
                writer=csv.writer(h,delimiter='\t');writer.writerow(['local_family','original_native_id'])
                for number,genes in groups(cluster_files[0]):
                    mapped=[ids[g] for g in genes]
                    if not mapped or len(set(mapped))!=len(mapped) or seen.intersection(mapped):raise ValueError('Non-unique family partition')
                    seen.update(mapped);family_count+=1;singletons+=len(mapped)==1
                    for original in sorted(mapped):writer.writerow([f'OG{number:07d}',original])
            if seen!=expected:raise ValueError('Native partition does not cover all input proteins')
            if list(native.glob('Trees_ids/*.txt')):raise ValueError('Unexpected gene-tree inference')
            for item in case['input_files']:
                if sha(inputs/item['filename'])!=item['sha256']:raise ValueError('Input changed during inference')
            artifacts={str(p.relative_to(folder)):sha(p) for p in [export,cluster_files[0],ids_file]}
            r=dict(status='complete_identity_checked_clade_discovery',clade_id=key,taxa=case['taxa'],proteins=len(expected),families=family_count,singleton_families=singletons,steps=steps,self_search=self_search,plan_sha256=plan_hash,artifacts=artifacts,interpretation='Native sequence-based family discovery within one guide-defined clade; exact protein coverage and identity verified. Not reconciled orthology, duplication inference or biological novelty.')
            save(receipt_path,r);completed.append(key);state.update(completed_clades=list(completed));save(root/'state.json',state)
            print(key,len(expected),family_count,'complete',flush=True)
        state.update(status='complete',completed_clades=completed,elapsed_seconds=time.time()-started);save(root/'state.json',state)
        save(root/'receipt.json',dict(status='complete_all_planned_clade_discovery',clades=len(completed),plan_sha256=plan_hash,clade_receipts={key:sha(root/key/'receipt.json') for key in completed},interpretation='Per-clade partitions completed; independent readback, guide-level singleton replacement/merging, gene trees and reconciliation remain required.'))
    except BaseException as exc:
        state.update(status='failed_requires_review',error=repr(exc),completed_clades=completed);save(root/'state.json',state);raise


if __name__=='__main__':main()
