#!/usr/bin/env python3
"""Install the official pinned Ubuntu archive locally and verify its release digest."""
import json, shutil, subprocess, tarfile, urllib.request
from pathlib import Path
from prepare_case_ancestral_neighborhoods import sha


def main():
    pp=Path('metadata/baliphy_install_plan_20260927.json');plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():assert sha(p)==h,p
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    archive=out/'bali-phy-4.3-ubuntu-24.04.tar.gz'
    with urllib.request.urlopen(plan['archive_url'],timeout=120) as response, archive.open('wb') as handle:shutil.copyfileobj(response,handle)
    assert sha(archive)==plan['archive_sha256']
    with tarfile.open(archive) as tar:tar.extractall(out/'install',filter='data')
    binaries=list((out/'install').glob('*/bin/bali-phy'));assert len(binaries)==1
    binary=binaries[0].resolve()
    for flag,name in [('--version','version.txt'),('--help','help.txt')]:
        result=subprocess.run([str(binary),flag],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=60)
        (out/name).write_text(result.stdout);assert result.returncode==0,(flag,result.returncode)
    assert '4.3' in (out/'version.txt').read_text()
    receipt=dict(status='official_archive_verified_binary_help_available_inference_untested',plan_sha256=sha(pp),archive_sha256=sha(archive),binary=str(binary),binary_sha256=sha(binary),artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},scope='Local user-owned install only. No system changes, project MCMC, capacity claim or convergence claim.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(str(binary))
if __name__=='__main__':main()
