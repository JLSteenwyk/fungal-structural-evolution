#!/usr/bin/env python3
"""Reproduce suspected exponent truncation using a pure installed native program."""
import json
from datetime import datetime,timezone
from pathlib import Path

from ancestral_chain_attempt import run_attempt,sha


def main():
    root=Path('data/software_audits/baliphy-native-double-format-cjson-20261003-v2').resolve()
    binary=Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/bali-phy').resolve()
    prlimit=Path('/usr/bin/prlimit');program=root/'format.hs'
    package=binary.parent.parent/'lib/bali-phy/haskell'
    paths=[binary,prlimit,program,*[package/name for name in ['Data/JSON.hs','Data/JSON/Encoding.hs','Data/JSON/Types/ToJSON.hs','Data/Text/Display.hs','Data/Text.hs']]]
    config=dict(command=[str(prlimit),'--as='+str(8*2**30),'--cpu=20','--fsize='+str(16*2**20),'--',str(binary),'run',str(program)],timeout_seconds=30,pins={str(p):sha(p) for p in paths})
    receipt=run_attempt(root/'native',config);result=json.loads(receipt.read_text());assert result['exit_code']==0
    output=receipt.parent/'stdout.log';lines=output.read_text().splitlines(); assert len(lines)==2
    values=json.loads(lines[0]); native_roundtrips=json.loads(lines[1])
    expected=[2.34e-10,2.34e-11,2.34e-20,2.34e-50,2.34e-100,2.34e-101,2.34e10,2.34e20,0.,1.,.25,4.]
    assert len(values)==len(expected)
    changed=[dict(index=i,expected=e,encoded=v,ratio=v/e) for i,(e,v) in enumerate(zip(expected,values)) if v!=e]
    import math
    assert all(math.isclose(v,e,rel_tol=2e-15,abs_tol=0) for v,e in zip(values,expected))
    assert native_roundtrips==[True]*12
    bindings=dict(config['pins']);bindings[str(Path(__file__))]=sha(__file__);bindings[str(receipt)]=sha(receipt)
    for name,digest in result['artifacts'].items():bindings[str(receipt.parent/name)]=digest
    proof=dict(status='passed_installed_native_cjson_fixed_double_roundtrip',checked_utc=datetime.now(timezone.utc).isoformat(),
        fixed_numeric_constants=12,native_exact_roundtrips=12,literal_parser_comparison_differences=changed,native_receipt=str(receipt),native_exit_code=0,
        source_hashes=bindings,scientific_eligibility=False,new_mcmc_runs=0,
        scope='Actual pure installed native CJSON encoding exactly roundtrips all12native fixed double values through native encode/read; cross-language literal differences remain explicit, including scientific exponents ending in zero. No installed source/binary mutation, native sampler restart or assertion relaxation. Full future model/logger transformation, all-role admission, native saved-frame validation and historical assessment remain required.')
    path=Path('metadata/baliphy_native_double_format_cjson_probe_20261003_v2.json')
    with path.open('x') as f:f.write(json.dumps(proof,indent=2)+'\n')
    print(json.dumps({k:v for k,v in proof.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
