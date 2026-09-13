#!/usr/bin/env python3
"""Acquire a pinned Pfam release and partition all profile HMMs for search."""
import datetime
import fcntl
import gzip
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = 'https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/'
EXPECTED = {'Pfam-A.hmm.gz': '7ab3c4e215d0daaea3004e37c4e24f8a',
            'Pfam-A.hmm.dat.gz': 'aab364689372b1846e8e880976e7f750',
            'Pfam-A.clans.tsv.gz': 'a70eeffcd8105ff094ea5834456d8bab',
            'Pfam.version.gz': '6c7beb9426eb0b1075e778c20b842b24',
            'active_site.dat.gz': 'f4161ba78f969238177654b32dbee4a2'}


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    root = ROOT / 'data/pfam/38.2'
    root.mkdir(parents=True, exist_ok=True)
    lock = (root / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    with urllib.request.urlopen(BASE + 'md5_checksums', timeout=90) as response:
        checksum_data = response.read()
    published = dict((line.split()[1], line.split()[0]) for line in checksum_data.decode().splitlines())
    if any(published.get(name) != value for name, value in EXPECTED.items()):
        raise ValueError('Publisher checksums differ from the pinned release')
    (root / 'md5_checksums').write_bytes(checksum_data)
    files = []
    for name, md5 in EXPECTED.items():
        path = root / name
        if not path.exists():
            temp = path.with_suffix('.partial')
            with urllib.request.urlopen(BASE + name, timeout=120) as response, temp.open('wb') as out:
                shutil.copyfileobj(response, out, length=1024 * 1024)
            if digest(temp, 'md5') != md5:
                raise ValueError('Downloaded file differs from publisher checksum')
            temp.replace(path)
        if digest(path, 'md5') != md5:
            raise ValueError('Changed cached Pfam input')
        target = path.with_suffix('')
        temp = target.with_suffix('.decompressing')
        with gzip.open(path, 'rb') as source, temp.open('wb') as out:
            shutil.copyfileobj(source, out, length=1024 * 1024)
        temp.replace(target)
        files.append({'file': name, 'url': BASE + name, 'md5': md5, 'sha256': digest(path),
                      'compressed_bytes': path.stat().st_size, 'uncompressed_file': target.name,
                      'uncompressed_sha256': digest(target), 'uncompressed_bytes': target.stat().st_size})
        print(name, 'verified', flush=True)
    version = (root / 'Pfam.version').read_text()
    if '38.2' not in version or '30134' not in version:
        raise ValueError('Unexpected Pfam version or family count')
    chunks = root / 'hmm_chunks'
    chunks.mkdir(exist_ok=True)
    handles = [(chunks / f'chunk-{i:03}.hmm.partial').open('w') for i in range(64)]
    names, lengths, counts = set(), [0] * 64, [0] * 64
    model, name, length, ga = [], None, None, False
    try:
        with (root / 'Pfam-A.hmm').open() as handle:
            for line in handle:
                model.append(line)
                if line.startswith('ACC '):
                    name = line.split()[1]
                elif line.startswith('LENG '):
                    length = int(line.split()[1])
                elif line.startswith('GA '):
                    ga = True
                elif line.strip() == '//':
                    if name is None or name in names or length is None or not ga:
                        raise ValueError('Invalid profile identity, length or gathering threshold')
                    names.add(name)
                    index = min(range(64), key=lambda i: (lengths[i], i))
                    handles[index].writelines(model)
                    lengths[index] += length
                    counts[index] += 1
                    model, name, length, ga = [], None, None, False
        if model or len(names) != 30134:
            raise ValueError('Truncated profile library or family count mismatch')
    finally:
        for handle in handles:
            handle.close()
    chunks_manifest = []
    for i, handle in enumerate(handles):
        temp = Path(handle.name)
        target = temp.with_suffix('')
        temp.replace(target)
        chunks_manifest.append({'path': str(target.relative_to(ROOT)), 'sha256': digest(target),
                                'profiles': counts[i], 'profile_match_states': lengths[i]})
    receipt = {'release': '38.2', 'families': len(names), 'version_text': version,
               'publisher_checksums_sha256': digest(root / 'md5_checksums'),
               'retrieved_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'files': files, 'chunks': chunks_manifest, 'script_sha256': digest(Path(__file__))}
    path = ROOT / 'metadata/pfam_release_receipt.json'
    if path.exists():
        previous = json.loads(path.read_text())
        if previous['files'] != files or previous['chunks'] != chunks_manifest:
            raise ValueError('Rebuilt Pfam library differs from previous receipt')
        receipt['retrieved_utc'] = previous['retrieved_utc']
    path.write_text(json.dumps(receipt, indent=2) + '\n')
    print('Verified and partitioned', len(names), 'profiles', flush=True)


if __name__ == '__main__':
    main()
