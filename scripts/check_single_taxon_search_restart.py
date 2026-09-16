#!/usr/bin/env python3
"""Test native one-taxon clustering from a complete synthetic self-search."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--fixture',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args(); out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=args.fixture/'2-taxon-no-fix/native-results'
    runs=list(source.glob('Results_*/WorkingDirectory'))
    if len(runs)!=1:raise ValueError('Ambiguous source')
    source=runs[0];wd=out/'WorkingDirectory';wd.mkdir()
    copied={}
    for name in ['Species0.fa','Blast0_0.txt.gz']:
        shutil.copyfile(source/name,wd/name);copied[name]=sha(wd/name)
    for name in ['SpeciesIDs.txt','SequenceIDs.txt']:
        prefix='0:' if name=='SpeciesIDs.txt' else '0_'
        (wd/name).write_text(''.join(line for line in (source/name).read_text().splitlines(keepends=True) if line.startswith(prefix)))
        copied[name]=sha(wd/name)
    root=Path(__file__).resolve().parents[1]
    command=[str(root/'.cache/envs/orthofinder/bin/orthofinder'),'-b',str(wd),'-t','2','-a','1','-S','diamond','--scores-v2','--only-groups','--no-fix-files','-I','1.2','-n','single_taxon_fixture']
    started=time.time();timeout=False
    with (out/'stdout.log').open('w') as log:
        p=subprocess.Popen(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'))
        try:rc=p.wait(timeout=180)
        except subprocess.TimeoutExpired:
            timeout=True;os.killpg(p.pid,signal.SIGTERM)
            try:rc=p.wait(timeout=10)
            except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);rc=p.wait()
    result=dict(status='observed_single_taxon_native_search_restart',returncode=rc,timed_out=timeout,elapsed_seconds=time.time()-started,command=command,
                original_input_sha256=copied,original_inputs_unchanged=all(sha(wd/n)==h for n,h in copied.items()),
                result_directories=[str(p.relative_to(out)) for p in out.rglob('Results_*') if p.is_dir()],
                fixture_receipt_sha256=sha(args.fixture/'receipt.json'),log_sha256=sha(out/'stdout.log'),script_sha256=sha(Path(__file__)),
                interpretation='Synthetic self-search copied from the two-taxon software fixture. Native restart is tested without adding a fictitious species. Exit status requires artifact validation; no biological inference or full-scale search launched.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
