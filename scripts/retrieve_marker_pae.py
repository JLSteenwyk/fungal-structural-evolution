#!/usr/bin/env python3
"""Retrieve version-matched PAE matrices for every model in a mapping snapshot."""
import argparse
import concurrent.futures
import datetime
import fcntl
import gzip
import hashlib
import json
import time
import urllib.request
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_pae(data, length):
    payload = json.loads(data)
    if not isinstance(payload, list) or len(payload) != 1:
        raise ValueError('Expected one monomer PAE object')
    matrix = np.asarray(payload[0]['predicted_aligned_error'], dtype=float)
    if matrix.shape != (length, length):
        raise ValueError(f'PAE dimensions {matrix.shape} do not match protein length {length}')
    if not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError('PAE must be finite and nonnegative')
    maximum = float(payload[0]['max_predicted_aligned_error'])
    # Some exports round matrix entries; permit half a unit of rounding.
    if not np.isfinite(maximum) or maximum < 0 or matrix.max() > maximum + .51:
        raise ValueError('PAE exceeds declared maximum')
    return matrix


def retrieve(model, cache):
    key = f"{model['model_id']}-v{model['version']}"
    url = model['pae_url']
    expected = f"https://alphafold.ebi.ac.uk/files/{model['model_id']}-predicted_aligned_error_v{model['version']}.json"
    if url != expected:
        raise ValueError('PAE URL does not match recorded model and version')
    path, receipt_path = cache / f'{key}.json.gz', cache / f'{key}.receipt.json'
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text())
        for field in ['model_id', 'version', 'length', 'sequence_sha256']:
            if receipt[field] != model[field]:
                raise ValueError('Cached PAE provenance mismatch')
        if receipt['url'] != url or sha(path) != receipt['gzip_sha256']:
            raise ValueError('Changed cached PAE')
        raw = gzip.decompress(path.read_bytes())
        if hashlib.sha256(raw).hexdigest() != receipt['json_sha256']:
            raise ValueError('Changed uncompressed PAE')
        validate_pae(raw, model['length'])
        return receipt
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'fungal-structural-evolution/1.0'})
            with urllib.request.urlopen(request, timeout=90) as response:
                raw = response.read()
            validate_pae(raw, model['length'])
            break
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    compressed = gzip.compress(raw, mtime=0)
    temporary = path.with_suffix('.tmp')
    temporary.write_bytes(compressed)
    temporary.replace(path)
    receipt = {field: model[field] for field in ['model_id', 'version', 'length', 'sequence_sha256']}
    receipt.update(status='verified', url=url, path=str(path.relative_to(ROOT)),
                   gzip_sha256=sha(path), json_sha256=hashlib.sha256(raw).hexdigest(),
                   compressed_bytes=len(compressed), json_bytes=len(raw),
                   retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    temporary = receipt_path.with_suffix('.tmp')
    temporary.write_text(json.dumps(receipt, indent=2) + '\n')
    temporary.replace(receipt_path)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = args.snapshot / 'receipt.json'
    mapping_receipt = json.loads(source.read_text())
    models_path = args.snapshot / 'model_provenance.json'
    if sha(models_path) != mapping_receipt['artifacts'][models_path.name]:
        raise ValueError('Changed source model provenance')
    models = json.loads(models_path.read_text())
    if len({(m['model_id'], m['version']) for m in models}) != len(models):
        raise ValueError('Duplicate model/version')
    if args.output.exists():
        raise FileExistsError('Choose a new immutable output snapshot')
    cache = ROOT / 'data/structures/pae'
    cache.mkdir(parents=True, exist_ok=True)
    lock = (cache / '.retrieval.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    args.output.mkdir(parents=True)
    records = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        futures = {executor.submit(retrieve, m, cache): m for m in models}
        for future in concurrent.futures.as_completed(futures):
            model = futures[future]
            try:
                record = future.result()
            except Exception as error:
                record = {'model_id': model['model_id'], 'version': model['version'],
                          'status': 'failed', 'error': str(error)}
            records.append(record)
            print(len(records), model['model_id'], record['status'], flush=True)
    records.sort(key=lambda r: (r['model_id'], r['version']))
    manifest = args.output / 'pae_manifest.json'
    manifest.write_text(json.dumps(records, indent=2) + '\n')
    verified = [r for r in records if r['status'] == 'verified']
    receipt = {'mapping_snapshot': str(args.snapshot), 'mapping_receipt_sha256': sha(source),
               'models_requested': len(models), 'models_verified': len(verified),
               'models_failed': len(records) - len(verified),
               'compressed_bytes': sum(r['compressed_bytes'] for r in verified),
               'json_bytes': sum(r['json_bytes'] for r in verified),
               'script_sha256': sha(Path(__file__)),
               'artifacts': {'pae_manifest.json': sha(manifest)}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == '__main__':
    main()
