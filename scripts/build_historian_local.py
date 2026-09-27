#!/usr/bin/env python3
"""Build pinned Historian with locally unpacked, versioned GSL dependencies."""
import json,subprocess
from pathlib import Path
from prepare_case_ancestral_neighborhoods import sha


def main():
    pp=Path('metadata/historian_build_plan_20260927.json');plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():assert sha(p)==h,p
    source=Path(plan['source']).resolve()
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip()==plan['commit']
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=source,text=True).strip()
    out=Path(plan['output']).resolve();out.mkdir(parents=True,exist_ok=False)
    deps=out/'packages';deps.mkdir();prefix=out/'dependencies';prefix.mkdir()
    with (out/'download.log').open('w') as h:
        subprocess.run(['apt-get','download']+plan['packages'],cwd=deps,stdout=h,stderr=subprocess.STDOUT,check=True)
    packages=sorted(deps.glob('*.deb'));assert len(packages)==3
    for p in packages:subprocess.run(['dpkg-deb','-x',str(p),str(prefix)],check=True)
    lib=prefix/'usr/lib/x86_64-linux-gnu'
    command=['make','-j4','CPP=g++','USING_BOOST=',f'GSL_FLAGS=-I{prefix}/usr/include',f'GSL_LIBS=-L{lib} -Wl,-rpath,{lib} -lgsl -lgslcblas','all']
    (out/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    (out/'compiler.txt').write_text(subprocess.check_output(['g++','--version'],text=True))
    with (out/'build.log').open('w') as h:subprocess.run(command,cwd=source,stdout=h,stderr=subprocess.STDOUT,check=True)
    binary=source/'bin/historian';assert binary.is_file()
    help_result=subprocess.run([str(binary),'-help'],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/'help.txt').write_text(help_result.stdout)
    assert 'recon' in help_result.stdout.lower()
    result=dict(status='historian_built_pending_functional_and_capacity_validation',commit=plan['commit'],plan_sha256=sha(pp),binary_path=str(binary),binary_sha256=sha(binary),help_exit_code=help_result.returncode,packages={p.name:sha(p) for p in packages},artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},scope='Local build only, no project inference or full-sampling capacity claim. GSL dependencies unpacked locally; no system packages installed. Source unchanged; compiler/flags recorded.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(result['status'],flush=True)

if __name__=='__main__':main()
