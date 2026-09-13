#!/usr/bin/env python3
"""Install the pinned official release and bundled tools in an isolated environment."""
import hashlib
import json
import tarfile
import subprocess
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
URL = 'https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.5/orthofinder-linux-intel-3.1.5.tar.gz'
EXPECTED = 'd74b5dbf9348e7ffdb735f1fda4ff9c49f504261811355519cf9dafcf80e7739'


def main():
    folder = ROOT / '.cache/software/orthofinder-3.1.5'
    archive = ROOT / '.cache/downloads/orthofinder-linux-intel-3.1.5.tar.gz'
    archive.parent.mkdir(parents=True, exist_ok=True)
    if not archive.exists():
        temp = archive.with_suffix('.partial')
        with urlopen(URL, timeout=90) as response, temp.open('wb') as out:
            for block in iter(lambda: response.read(1024 * 1024), b''):
                out.write(block)
        if hashlib.sha256(temp.read_bytes()).hexdigest() != EXPECTED:
            raise ValueError('Release checksum mismatch')
        temp.replace(archive)
    if hashlib.sha256(archive.read_bytes()).hexdigest() != EXPECTED:
        raise ValueError('Cached release checksum mismatch')
    if not folder.exists():
        folder.mkdir(parents=True)
        with tarfile.open(archive) as handle:
            handle.extractall(folder, filter='data')
    packages = list(folder.glob('*/pyproject.toml'))
    if len(packages) != 1:
        raise ValueError('Expected one package root')
    prefix = ROOT / '.cache/envs/orthofinder'
    if not (prefix / 'bin/python').exists():
        subprocess.run(['conda', 'create', '--yes', '--solver', 'libmamba', '--prefix', str(prefix),
                        '--override-channels', '-c', 'conda-forge', 'python=3.12', 'pip'], check=True)
    dependency_lock = ROOT / 'environments/orthofinder-pip.lock.txt'
    if dependency_lock.exists():
        subprocess.run([str(prefix / 'bin/python'), '-m', 'pip', 'install', '-r', str(dependency_lock)], check=True)
    subprocess.run([str(prefix / 'bin/python'), '-m', 'pip', 'install', str(packages[0].parent)], check=True)
    freeze = subprocess.run([str(prefix / 'bin/python'), '-m', 'pip', 'freeze'], capture_output=True, text=True, check=True).stdout
    dependency_lock.write_text('\n'.join(line for line in freeze.splitlines() if not line.lower().startswith('orthofinder')) + '\n')
    executable = prefix / 'bin/orthofinder'
    (ROOT / 'metadata/orthofinder_software_receipt.json').write_text(json.dumps({
        'version': '3.1.5', 'url': URL, 'archive_sha256': EXPECTED,
        'checksum_source': 'GitHub official release asset API digest',
        'executable_path': str(executable.relative_to(ROOT)),
        'executable_sha256': hashlib.sha256(executable.read_bytes()).hexdigest()}, indent=2) + '\n')
    print(executable)


if __name__ == '__main__':
    main()
