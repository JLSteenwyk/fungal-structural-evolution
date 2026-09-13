#!/usr/bin/env python3
"""Pin official ESMFold v1 files and verify cached weights against publisher SHA256."""
import argparse
import hashlib
import json
import os
import urllib.request
from pathlib import Path
from prepare_pfam import ROOT, digest

REVISION = 'ba837a39b67e59941c3f017d6c2a064f567038d9'
MODEL = 'facebook/esmfold_v1'
WEIGHT_SHA = '9a865162cdcaac8d5385c908fd0c620fccf0405adb3a5b49566a8d88831977ac'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cached-weights', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable checkpoint directory')
    api = f'https://huggingface.co/api/models/{MODEL}/revision/{REVISION}?blobs=true'
    with urllib.request.urlopen(api, timeout=90) as response:
        raw = response.read()
    info = json.loads(raw)
    entries = {r['rfilename']: r for r in info['siblings']}
    if info['sha'] != REVISION or entries['model.safetensors']['lfs']['sha256'] != WEIGHT_SHA:
        raise ValueError('Publisher checkpoint identity changed')
    if digest(args.cached_weights) != WEIGHT_SHA:
        raise ValueError('Cached weights do not match publisher checksum')
    args.output.mkdir(parents=True)
    (args.output / 'publisher.json').write_bytes(raw)
    sources = {}
    for name in ['config.json', 'vocab.txt', 'special_tokens_map.json', 'tokenizer_config.json']:
        url = f'https://huggingface.co/{MODEL}/resolve/{REVISION}/{name}'
        with urllib.request.urlopen(url, timeout=90) as response:
            data = response.read()
        blob = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        if len(data) != entries[name]['size'] or blob != entries[name]['blobId']:
            raise ValueError('Publisher Git blob identity mismatch: ' + name)
        (args.output / name).write_bytes(data)
        sources[name] = url
    # Reuse the verified local immutable blob; no second 8.4 GB copy is needed.
    os.symlink(args.cached_weights.resolve(), args.output / 'model.safetensors')
    receipt = {'model': MODEL, 'revision': REVISION, 'publisher_api': api,
               'weights_sha256': WEIGHT_SHA, 'weights_bytes': args.cached_weights.stat().st_size,
               'weights_url': f'https://huggingface.co/{MODEL}/resolve/{REVISION}/model.safetensors',
               'sources': sources, 'script_sha256': digest(Path(__file__)),
               'artifacts': {p.name: digest(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
