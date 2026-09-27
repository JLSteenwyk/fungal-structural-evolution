#!/usr/bin/env python3
"""Run upstream Historian tests with the audited local build configuration."""
import json, subprocess, time
from pathlib import Path
from prepare_case_ancestral_neighborhoods import sha

def main():
    pp=Path('metadata/historian_test_plan_20260927.json')
    plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items(): assert sha(p)==h,p
    build=Path(plan['build'])
    source=Path(plan['source']).resolve()
    receipt=json.loads((build/'receipt.json').read_text())
    assert sha(receipt['binary_path'])==receipt['binary_sha256']
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    command=json.loads((build/'command.json').read_text())[:-1]+['-k','test']
    (out/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    start=time.time()
    with (out/'tests.log').open('w') as h:
        try:
            result=subprocess.run(command,cwd=source,stdout=h,stderr=subprocess.STDOUT,timeout=plan['timeout_seconds'])
            code=result.returncode;status='upstream_tests_passed' if code==0 else 'upstream_test_failures_require_review'
        except subprocess.TimeoutExpired:
            code=None;status='upstream_tests_timeout_require_review'
    assert sha(receipt['binary_path'])==receipt['binary_sha256']
    result=dict(status=status,exit_code=code,elapsed_seconds=time.time()-start,plan_sha256=sha(pp),binary_sha256=receipt['binary_sha256'],artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},scope='Upstream test suite only; does not validate project root placement, full-family capacity or posterior uncertainty.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(status,flush=True)
if __name__=='__main__': main()
