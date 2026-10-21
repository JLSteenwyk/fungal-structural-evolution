#!/usr/bin/env python3
"""Retry only explicit full-atlas PAE retrieval failures after a completed recovery."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import time

from retrieve_marker_pae import retrieve
from prepare_full_atlas_pae_failure_retry_v1 import prepare


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def key(model):
    return f"{model['model_id']}-v{model['version']}"


def submit_window(models, executor, size):
    models = iter(models); pending = {}; exhausted = False
    while pending or not exhausted:
        while not exhausted and len(pending) < size:
            try: model = next(models)
            except StopIteration: exhausted = True; break
            pending[executor.submit(retrieve, model)] = model
        if not pending: break
        done, _ = wait(pending, return_when='FIRST_COMPLETED')
        for future in done: yield pending.pop(future), future


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--queue', type=Path, required=True)
    p.add_argument('--manifest', type=Path, action='append', required=True)
    p.add_argument('--cache', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    p.add_argument('--http-workers', type=int, default=1)
    p.add_argument('--future-window', type=int, default=1)
    p.add_argument('--emergency-free-gib', type=int, default=512)
    p.add_argument('--max-new-cache-bytes', type=int, default=1099511627776)
    a = p.parse_args()
    if a.output.exists() or a.receipt.exists(): raise FileExistsError('Fresh output and receipt paths required')
    if a.http_workers != 1 or a.future_window != 1: raise ValueError('Retry stage is intentionally serial')
    if not a.cache.is_dir(): raise FileNotFoundError('Preserved PAE cache is required')
    if shutil.disk_usage(a.cache).free < a.emergency_free_gib * 2**30: raise RuntimeError('Insufficient initial disk reserve')
    a.output.mkdir(parents=True)
    preparation = prepare(a.queue, a.manifest, a.output / 'queue-preparation', a.output / 'queue-preparation.json', Path(__file__).with_name('prepare_full_atlas_pae_failure_retry_v1.py'))
    retry_queue = Path(preparation['retry_queue'])
    counts = Counter(attempted=0, verified=0, failed=0, compressed_bytes=0, json_bytes=0)
    errors = Counter(); started = time.monotonic(); manifest = a.output / 'retry_dispositions.jsonl.gz'

    def state(stage):
        value = dict(stage=stage, counts=dict(counts), error_types=dict(errors), expected_models=preparation['retry_models'],
                     elapsed_seconds=time.monotonic()-started, http_workers=1, future_window=1)
        tmp = a.output / 'state.tmp'; tmp.write_text(json.dumps(value, indent=2)+'\n'); tmp.replace(a.output/'state.json')

    def fetch(model):
        if shutil.disk_usage(a.cache).free < a.emergency_free_gib * 2**30: raise RuntimeError('Emergency disk reserve reached')
        return retrieve(model, a.cache / model['sequence_sha256'][:2])

    state('retrying_explicit_failures')
    with manifest.open('xb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as archive:
        with io.TextIOWrapper(archive, encoding='utf-8', newline='') as out, retry_queue.open() as handle, ThreadPoolExecutor(max_workers=1) as executor:
            models = (json.loads(line) for line in handle)
            for model, future in submit_window(models, executor, 1):
                counts['attempted'] += 1
                try: result = future.result()
                except Exception as error:
                    if isinstance(error, RuntimeError) and str(error) == 'Emergency disk reserve reached': raise
                    counts['failed'] += 1; errors[type(error).__name__] += 1
                    row = dict(model=model, status='retrieval_failed', error_type=type(error).__name__, error=str(error))
                else:
                    if result.get('status') != 'verified' or any(result.get(f) != model[f] for f in ('model_id','version','sequence_sha256','length')):
                        raise ValueError('Retrieved receipt differs from retry target')
                    receipt_path = Path(result['path']).with_name(f'{key(model)}.receipt.json')
                    if json.loads(receipt_path.read_text()) != result: raise ValueError('Retained receipt differs from retriever result')
                    counts['verified'] += 1; counts['compressed_bytes'] += result['compressed_bytes']; counts['json_bytes'] += result['json_bytes']
                    row = dict(model=model, status='retrieval_verified', receipt_path=str(receipt_path), receipt_sha256=sha(receipt_path), receipt=result)
                out.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
                if counts['attempted'] % 25 == 0: out.flush(); state('retrying_explicit_failures')
                if counts['compressed_bytes'] > a.max_new_cache_bytes: raise RuntimeError('Declared additional-cache allowance reached')
    if counts['attempted'] != preparation['retry_models'] or counts['attempted'] != counts['verified'] + counts['failed']:
        raise ValueError('Retry queue census differs')
    result = dict(status='completed_full_atlas_pae_explicit_failure_retry_pending_independent_readback',
                  checked_utc=datetime.now(timezone.utc).isoformat(), counts=dict(counts), error_types=dict(errors),
                  expected_models=preparation['retry_models'], queue_preparation=str(a.output/'queue-preparation.json'),
                  queue_preparation_sha256=sha(a.output/'queue-preparation.json'), retry_manifest=str(manifest), retry_manifest_sha256=sha(manifest),
                  source_hashes={str(p):sha(p) for p in [Path(__file__), Path(__file__).with_name('prepare_full_atlas_pae_failure_retry_v1.py'), Path(__file__).with_name('retrieve_marker_pae.py'), a.queue,*a.manifest]}, scientific_eligibility=False,
                  gpu=False, new_predictions=0,
                  scope='A serial retry of only explicit retained failures after a completed source-bound recovery. It preserves cache provenance and requires a separate combined independent reader.')
    a.receipt.write_text(json.dumps(result, indent=2)+'\n')
    state('all_explicit_failure_retries_complete_pending_independent_readback')
    print(json.dumps({k:v for k,v in result.items() if k != 'source_hashes'}, indent=2))

if __name__ == '__main__': main()
