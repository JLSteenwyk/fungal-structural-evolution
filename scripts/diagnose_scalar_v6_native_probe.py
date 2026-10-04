#!/usr/bin/env python3
"""Bounded pure-native isolation and GDB inspection of the failed numeric fixture."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

from ancestral_chain_attempt import run_attempt,sha
from reference_measurement_union_sources import verify


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);a=p.parse_args();assert not a.receipt.exists()
    root=a.output.resolve();root.mkdir(exist_ok=False)
    original=Path('data/software_audits/baliphy-scalar-json-v6-20261004-v5/value-probe.hs').resolve()
    text=original.read_text();head=text.split('main = do\n')[0]
    datafile=Path('data/software_audits/baliphy-scalar-json-v6-20261004-v5/nonfinite-values.txt').resolve()
    binary=Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/bali-phy').resolve()
    sources=[original,datafile,Path(__file__).resolve(),binary,Path('/usr/bin/prlimit'),Path('/usr/bin/gdb')]
    libraries=list((binary.parent.parent/'lib/bali-phy/haskell').rglob('*.hs'))
    pins={str(q):sha(q) for q in sources+libraries};cases=[]
    variants={
      'gdb-original':text,
      'ordinary-quality':head+'''main = do
  let simple = [(J.toJSONKey "finite", J.Array (map J.FNumber (take 12 probeValuesV6)))]
  T.putStrLn (J.fromEncoding (projectV6EncodeRecord 0 (projectV6ContextValue simple) simple))
''',
      'runtime-nonfinite':head+'''main = do
  contents <- T.readFile '''+json.dumps(str(datafile))+'''
  let special = map read (lines (Text.unpack contents)) :: [Double]
  T.putStrLn (J.cjsonToText (J.toCJSON (map isInfinite special)))
  T.putStrLn (J.cjsonToText (J.toCJSON (map isNaN special)))
  T.putStrLn (J.cjsonToText (J.toCJSON special))
''',
      'extreme-finite':head+'''main = do
  T.putStrLn (J.cjsonToText (J.toCJSON (drop 15 probeValuesV6)))
'''}
    for name,source in variants.items():
        program=root/(name+'.hs');program.write_text(source)
        command=['/usr/bin/prlimit','--as='+str(12*2**30),'--cpu=120','--fsize='+str(16*2**20),'--']
        if name=='gdb-original':command+=['/usr/bin/gdb','--batch','--quiet','-ex','run','-ex','bt 25','--args']
        command += [str(binary),'run',str(program)]
        config=dict(command=command,timeout_seconds=180,pins={**pins,str(program):sha(program)})
        target=root/'native'/name;assert not target.exists()
        receipt=run_attempt(target,config);r=json.loads(receipt.read_text())
        pins.update(config['pins']);pins[str(receipt)]=sha(receipt)
        pins[str(target/'configuration.json')]=sha(target/'configuration.json')
        for q,d in r['artifacts'].items():pins[str(receipt.parent/q)]=d
        cases.append(dict(case=name,exit_code=r['exit_code'],status=r['status'],native_receipt=str(receipt)))
        print(json.dumps(cases[-1]),flush=True)
    verify(pins)
    result=dict(status='completed_scalar_v6_native_fixture_isolation_not_qualification',
        checked_utc=datetime.now(timezone.utc).isoformat(),cases=cases,source_hashes=pins,
        native_mcmc_runs=0,software_fix_qualified=False,installed_sources_changed=False,
        original_failed_attempts_restarted=False,gpu=False,scientific_eligibility=False,
        scope='New pure-native and debugger programs isolate the original numeric fixture failure. '
              'GDB exit is not native success. Original failed program and output bytes remain unchanged. '
              'No successful full logger qualification, production sampler, crash repair or biological acceptance.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')


if __name__=='__main__':main()
