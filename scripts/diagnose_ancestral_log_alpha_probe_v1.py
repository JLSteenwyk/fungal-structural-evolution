#!/usr/bin/env python3
"""Inspect the failed pure-native probe under an unchanged child resource cap."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import run_attempt, sha
from reference_measurement_union_sources import bind, verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args()
    assert not args.output.exists() and not args.receipt.exists()
    original=Path('results/ancestral/full-current-log-alpha-diagnostic-20261005-v1/deterministic-native-probe')
    config_path=original/'configuration.json'
    previous=json.loads(config_path.read_text())
    failed_path=original/'attempt-0001/receipt.json'
    failed=json.loads(failed_path.read_text())
    assert failed['exit_code']==-11 and failed['status']=='failed'
    pins=dict(previous['pins'])
    verify(pins)
    bind(pins,config_path);bind(pins,failed_path)
    for name,digest in failed['artifacts'].items():bind(pins,failed_path.parent/name,digest)
    bind(pins,Path('/usr/bin/gdb'));bind(pins,Path(__file__))
    command=previous['command']
    assert command[:5]==['/usr/bin/prlimit','--as='+str(8*2**30),'--cpu=60','--fsize='+str(16*2**20),'--']
    config=dict(command=['/usr/bin/prlimit','--as='+str(12*2**30),'--cpu=120',
        '--fsize='+str(128*2**20),'--','/usr/bin/gdb','--batch','--quiet','--nx',
        '-ex','run','-ex','bt 25','-ex','info registers rip rsp rbp','--args',*command],
        timeout_seconds=180,pins={str(Path(path).resolve()):digest for path,digest in pins.items()})
    attempt=run_attempt(args.output,config)
    result=json.loads(attempt.read_text())
    stdout=attempt.parent/'stdout.log'
    text=stdout.read_text()
    assert result['exit_code']==0 and 'SIGSEGV' in text
    for path,digest in config['pins'].items():bind(pins,path,digest)
    bind(pins,attempt)
    for name,digest in result['artifacts'].items():bind(pins,attempt.parent/name,digest)
    verify(pins)
    proof=dict(status='captured_original_deterministic_probe_native_signal',checked_utc=datetime.now(timezone.utc).isoformat(),
        original_native_exit_code=-11,gdb_exit_code=0,signal='SIGSEGV',source_hashes=pins,
        scientific_eligibility=False,new_mcmc_runs=0,original_failed_attempt_restarted=False,gpu=False,
        scope='Separate debugger inspection of the exact failed deterministic native program with '
              'original child limits preserved. Debugger exit zero is a captured signal, not native '
              'probe success, a chain replay, repair or posterior acceptance. Original failed bytes retained.')
    with args.receipt.open('x') as handle:json.dump(proof,handle,indent=2);handle.write('\n')
    print(json.dumps({k:v for k,v in proof.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
