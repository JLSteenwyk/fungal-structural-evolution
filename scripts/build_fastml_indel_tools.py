#!/usr/bin/env python3
"""Build the archived FastML indel tools with recorded source and compiler."""
import hashlib,json,subprocess,time
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    pp=Path('metadata/fastml_indel_build_plan_20260927.json');p=json.loads(pp.read_text())
    for name,h in p['pins'].items():assert sha(name)==h,name
    root=Path(p['source']).resolve();out=Path(p['output']);out.mkdir(parents=True,exist_ok=False)
    source={str(f.relative_to(root)):sha(f) for f in root.rglob('*') if f.is_file()};(out/'source_files.json').write_text(json.dumps(source,indent=2)+'\n')
    compiler=subprocess.check_output(['g++','--version'],text=True);(out/'compiler.txt').write_text(compiler)
    commands=[['make','-j4','CXX=g++ -std=gnu++11','libs'],['make','-j4','CXX=g++ -std=gnu++11','-C','programs','indelCoder','gainLoss']]
    with (out/'build.log').open('w') as h:
        for command in commands:subprocess.run(command,cwd=root,stdout=h,stderr=subprocess.STDOUT,check=True)
    binaries={}
    for name in ['indelCoder','gainLoss']:
        path=root/'programs'/name/name;assert path.is_file();binaries[name]=dict(path=str(path),sha256=sha(path))
    result=dict(status='built_fastml_311_indel_tools_pending_functional_validation',plan_sha256=sha(pp),commands=commands,binaries=binaries,source_manifest_sha256=sha(out/'source_files.json'),compiler_sha256=sha(out/'compiler.txt'),build_log_sha256=sha(out/'build.log'),scope='Build only. No indel inference, capacity claim or scientific qualification.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

if __name__=='__main__':main()
