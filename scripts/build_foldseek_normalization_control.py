#!/usr/bin/env python3
"""Build original and minimally patched pinned Foldseek in an isolated directory."""
import argparse
import hashlib
import json
import shutil
import subprocess
import tarfile
import urllib.request
from pathlib import Path

REVISION = 'e3fadcd07f971e864c094ac4f3a78bf4ed845e07'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--cmake', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new isolated output')
    a.output = a.output.resolve(); a.cmake = a.cmake.resolve()
    a.output.mkdir(parents=True)
    url = f'https://codeload.github.com/steineggerlab/foldseek/tar.gz/{REVISION}'
    archive = a.output / 'source.tar.gz'
    with urllib.request.urlopen(url, timeout=120) as response, archive.open('wb') as f:
        shutil.copyfileobj(response, f)
    with tarfile.open(archive) as tar:
        tar.extractall(a.output, filter='data')
    source = a.output / ('foldseek-' + REVISION)
    build = a.output / 'build'
    command = [str(a.cmake), '-S', str(source), '-B', str(build), '-DCMAKE_BUILD_TYPE=Release',
               '-DHAVE_AVX2=1', '-DENABLE_PROSTT5=0', '-DENABLE_CUDA=0']
    config = {'revision': REVISION, 'source_url': url, 'archive_sha256': sha(archive),
              'producer_sha256': sha(Path(__file__)), 'configure_command': command,
              'cmake_version': subprocess.check_output([str(a.cmake), '--version'], text=True),
              'compiler_version': subprocess.check_output(['g++', '--version'], text=True),
              'rust_version': subprocess.check_output(['rustc', '--version'], text=True),
              'interpretation': 'Same source and toolchain for unmodified and output-normalization patch builds; production binary unchanged. No biological result until rescoring readback.'}
    (a.output / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
    with (a.output / 'configure.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    compile_command = [str(a.cmake), '--build', str(build), '--target', 'foldseek', '--parallel', '4']
    with (a.output / 'unmodified-build.log').open('w') as log:
        subprocess.run(compile_command, stdout=log, stderr=subprocess.STDOUT, check=True)
    binary = build / 'src/foldseek'
    shutil.copy2(binary, a.output / 'foldseek-unmodified')
    target = source / 'src/strucclustutils/structureconvertalis.cpp'
    original = target.read_text()
    block = original.split('case LocalParameters::OUTFMT_ALNTMSCORE:')[1].split('break;')[0]
    old = 'std::min(res.qEndPos - res.qStartPos, res.dbEndPos - res.dbStartPos)'
    new = 'std::min(res.qEndPos - res.qStartPos + 1, res.dbEndPos - res.dbStartPos + 1)'
    if block.count(old) != 1 or original.count(block) != 1:
        raise ValueError('Expected unique output normalization block absent')
    updated = original.replace(block, block.replace(old, new))
    target.write_text(updated)
    import difflib
    (a.output / 'normalization.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True), updated.splitlines(True), fromfile='a/src/strucclustutils/structureconvertalis.cpp', tofile='b/src/strucclustutils/structureconvertalis.cpp')))
    with (a.output / 'patched-build.log').open('w') as log:
        subprocess.run(compile_command, stdout=log, stderr=subprocess.STDOUT, check=True)
    shutil.copy2(binary, a.output / 'foldseek-inclusive-span')
    receipt = {'status': 'complete_isolated_original_and_patched_build', 'config_sha256': sha(a.output / 'config.json'),
               'compile_command': compile_command, 'unmodified_source_sha256': hashlib.sha256(original.encode()).hexdigest(),
               'patched_source_sha256': sha(target),
               'artifacts': {name: sha(a.output / name) for name in ['foldseek-unmodified', 'foldseek-inclusive-span', 'normalization.patch']}}
    (a.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
