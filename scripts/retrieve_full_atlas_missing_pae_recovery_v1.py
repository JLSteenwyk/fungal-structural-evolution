#!/usr/bin/env python3
"""Retrieve the source-bound unfinished portion of an interrupted PAE queue."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
import gzip
import io
import json
from pathlib import Path
import shutil
import time

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify
from retrieve_marker_pae import retrieve


def completed_requests(models, executor, submit, window):
    models = iter(models)
    pending, exhausted = {}, False
    while pending or not exhausted:
        while not exhausted and len(pending) < window:
            try:
                model = next(models)
            except StopIteration:
                exhausted = True
                break
            pending[executor.submit(submit, model)] = model
        if not pending:
            break
        ready, _ = wait(pending, return_when=FIRST_COMPLETED)
        for future in ready:
            yield pending.pop(future), future


def model_key(model):
    return f"{model['model_id']}-v{model['version']}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--preflight', action='store_true', help='Validate immutable inputs without network or writes')
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError('Fresh receipt path required')
    plan = json.loads(args.plan.read_text())
    pins = dict(plan['pins'])
    verify(pins)
    preparation_path = Path(plan['recovery_preparation'])
    preparation = json.loads(preparation_path.read_text())
    if preparation.get('status') != 'prepared_interrupted_full_atlas_pae_recovery_queue':
        raise ValueError('Source-bound recovery preparation receipt required')
    queue = Path(preparation['recovery_queue'])
    if sha(queue) != preparation['recovery_queue_sha256']:
        raise ValueError('Recovery queue digest differs')
    if preparation['recovery_models'] != plan['expected_models']:
        raise ValueError('Plan model count differs from recovery preparation')
    root, cache = Path(plan['output']), Path(plan['cache']).resolve()
    if root.exists():
        raise FileExistsError('Fresh recovery output root required')
    if not cache.is_dir():
        raise FileNotFoundError('Preserved source cache is required')
    if shutil.disk_usage(cache).free < plan['resources']['emergency_free_disk_gib'] * 2**30:
        raise RuntimeError('Insufficient initial emergency disk reserve')
    if args.preflight:
        observed = 0
        with queue.open() as handle:
            for line in handle:
                model = json.loads(line)
                directory = cache / model['sequence_sha256'][:2]
                if (directory / f"{model_key(model)}.receipt.json").exists():
                    raise RuntimeError('Recovery target already has a retained receipt')
                observed += 1
        if observed != plan['expected_models']:
            raise ValueError('Recovery queue count differs from plan')
        print(json.dumps(dict(status='passed_full_atlas_missing_pae_recovery_preflight',
                              expected_models=observed, cache=str(cache),
                              free_bytes=shutil.disk_usage(cache).free,
                              scope='Read-only source, cache-collision and disk-reserve validation.'), indent=2))
        return
    root.mkdir(parents=True)
    counts = Counter(attempted=0, verified=0, failed=0, compressed_bytes=0, json_bytes=0)
    errors = Counter()
    started = time.monotonic()
    manifest = root / 'recovery_dispositions.jsonl.gz'

    def fetch(model):
        if shutil.disk_usage(cache).free < plan['resources']['emergency_free_disk_gib'] * 2**30:
            raise RuntimeError('Emergency disk reserve reached')
        directory = cache / model['sequence_sha256'][:2]
        directory.mkdir(exist_ok=True)
        expected_receipt = directory / f"{model_key(model)}.receipt.json"
        if expected_receipt.exists():
            raise RuntimeError('Recovery target already has a retained receipt')
        return retrieve(model, directory)

    def state(stage):
        value = dict(stage=stage, counts=dict(counts), error_types=dict(errors),
                     expected_models=plan['expected_models'], elapsed_seconds=time.monotonic() - started,
                     http_workers=plan['resources']['http_workers'], future_window=plan['future_window'])
        temporary = root / 'state.tmp'
        temporary.write_text(json.dumps(value, indent=2) + '\n')
        temporary.replace(root / 'state.json')

    state('retrieving_recovery_queue')
    with manifest.open('xb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as archive:
        with io.TextIOWrapper(archive, encoding='utf-8', newline='') as output, queue.open() as handle, ThreadPoolExecutor(max_workers=plan['resources']['http_workers']) as executor:
            models = (json.loads(line) for line in handle)
            seen = set()
            for model, future in completed_requests(models, executor, fetch, plan['future_window']):
                identifier = model_key(model)
                if identifier in seen:
                    raise ValueError('Duplicate recovery-queue model')
                seen.add(identifier)
                counts['attempted'] += 1
                try:
                    result = future.result()
                except Exception as error:
                    if isinstance(error, RuntimeError) and str(error) == 'Emergency disk reserve reached':
                        raise
                    errors[type(error).__name__] += 1
                    counts['failed'] += 1
                    row = dict(model=model, status='retrieval_failed', error_type=type(error).__name__, error=str(error))
                else:
                    if result.get('status') != 'verified' or any(result.get(field) != model[field] for field in ('model_id', 'version', 'sequence_sha256', 'length')):
                        raise ValueError('Retrieved receipt differs from recovery target')
                    receipt_path = Path(result['path']).with_name(f"{identifier}.receipt.json")
                    if json.loads(receipt_path.read_text()) != result:
                        raise ValueError('Retained receipt differs from retriever result')
                    counts['verified'] += 1
                    counts['compressed_bytes'] += result['compressed_bytes']
                    counts['json_bytes'] += result['json_bytes']
                    row = dict(model=model, status='retrieval_verified', receipt_path=str(receipt_path),
                               receipt_sha256=sha(receipt_path), receipt=result)
                output.write(json.dumps(row, separators=(',', ':'), allow_nan=False) + '\n')
                if counts['attempted'] % 100 == 0:
                    output.flush()
                    state('retrieving_recovery_queue')
                if counts['compressed_bytes'] > plan['resources']['maximum_new_cache_bytes']:
                    raise RuntimeError('Declared additional-cache allowance reached')
    if counts['attempted'] != plan['expected_models'] or counts['attempted'] != counts['verified'] + counts['failed']:
        raise ValueError('Recovery queue census differs from plan')
    verify(pins)
    for path in (args.plan, preparation_path, queue, manifest, Path(__file__)):
        bind(pins, path)
    verify(pins)
    result = dict(status='completed_full_atlas_missing_pae_recovery_attempts_pending_independent_readback',
                  checked_utc=datetime.now(timezone.utc).isoformat(), counts=dict(counts), error_types=dict(errors),
                  expected_models=plan['expected_models'], source_hashes=pins, scientific_eligibility=False,
                  gpu=False, new_predictions=0,
                  scope='Each previously failed or unattempted source-bound model was attempted once under the '
                        'qualified retriever. Existing verified cache entries are never rewritten. Successful and '
                        'failed recovery dispositions remain explicit; combined-cache acceptance requires a separate '
                        'order-independent independent reader.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    state('all_recovery_attempts_complete_pending_independent_readback')
    print(json.dumps({key: value for key, value in result.items() if key != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
