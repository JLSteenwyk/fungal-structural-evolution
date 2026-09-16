#!/usr/bin/env python3
"""Exercise native discovery on deterministic, non-biological software fixtures."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import signal
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); out=args.output.resolve(); out.mkdir(parents=True,exist_ok=False)
    rng=random.Random(20260916); alphabet='ACDEFGHIKLMNPQRSTVWY'
    families=[''.join(rng.choice(alphabet) for _ in range(100+15*g)) for g in range(24)]
    results=[]
    for taxa,no_fix in [(1,False),(2,False),(1,True),(2,True)]:
        case=out/(f'{taxa}-taxon-'+('no-fix' if no_fix else 'default')); source=case/'inputs'; source.mkdir(parents=True)
        expected=[]
        for sid in range(taxa):
            records=[]
            for family,sequence in enumerate(families):
                # Exact synthetic pairs test grouping and original-ID preservation.
                for copy in range(2):
                    name=f'{100+sid}_{family*2+copy}'
                    records.append(f'>{name}\n{sequence}\n'); expected.append(name)
            (source/f'SyntheticTaxon{sid}.faa').write_text(''.join(records))
        before={p.name:sha(p) for p in source.iterdir()}
        command=[str(ROOT/'.cache/envs/orthofinder/bin/orthofinder'),'-f',str(source),'-o',str(case/'native-results'),
                 '-t','2','-a','1','-S','diamond','--scores-v2','--only-groups','-I','1.2']
        if no_fix: command.append('--no-fix-files')
        started=time.time(); timed_out=False
        with (case/'stdout.log').open('w') as log:
            p=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
                               env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'))
            try: rc=p.wait(timeout=180)
            except subprocess.TimeoutExpired:
                timed_out=True; os.killpg(p.pid,signal.SIGTERM)
                try: rc=p.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid,signal.SIGKILL);rc=p.wait()
        after={p.name:sha(p) for p in source.iterdir() if p.is_file()}
        result=dict(case=case.name,taxa=taxa,no_fix_files=no_fix,expected_proteins=len(expected),expected_ids=expected,command=command,
                    returncode=rc,timed_out=timed_out,elapsed_seconds=time.time()-started,
                    source_unchanged=before==after,source_sha256=before,
                    result_directories=[str(p.relative_to(out)) for p in (case/'native-results').glob('Results_*')],
                    log_sha256=sha(case/'stdout.log'))
        results.append(result); print(taxa,'taxa',rc,'timeout',timed_out,flush=True)
    receipt=dict(status='completed_native_clade_discovery_software_observation',cases=results,
                 script_sha256=sha(Path(__file__)),executable_sha256=sha(ROOT/'.cache/envs/orthofinder/bin/orthofinder'),
                 interpretation='Synthetic random amino-acid strings and exact copies exercise software behavior and ID preservation. No real biological subset or pilot; no claim of homology/orthology accuracy, sensitivity or full-scale resource validation. Outputs require independent readback; exit zero alone is not success.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
