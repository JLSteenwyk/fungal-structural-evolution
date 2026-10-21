#!/usr/bin/env python3
"""Independently audit preserved and recovery PAE dispositions by source key.

Unlike the superseded reader, this reader does not assume that a concurrent
downloader completes records in source-queue order.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import tempfile


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def model_key(model):
    return f"{model['model_id']}-v{model['version']}"


def validate_matrix(raw, model):
    payload = json.loads(raw)
    if not isinstance(payload, list) or len(payload) != 1 or not isinstance(payload[0], dict):
        raise ValueError('Invalid PAE JSON container')
    body = payload[0]
    values, maximum = body.get('predicted_aligned_error'), body.get('max_predicted_aligned_error')
    if not isinstance(values, list) or len(values) != model['length'] or not isinstance(maximum, (int, float)):
        raise ValueError('Invalid PAE dimensions or declared maximum')
    if not math.isfinite(float(maximum)) or float(maximum) <= 0:
        raise ValueError('Invalid declared PAE maximum')
    entries = 0
    for row in values:
        if not isinstance(row, list) or len(row) != model['length']:
            raise ValueError('Non-square PAE matrix')
        for value in row:
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)) or value < 0:
                raise ValueError('Non-finite or negative PAE value')
            if value > maximum + .51:
                raise ValueError('PAE value exceeds declared export bound')
            entries += 1
    return entries


def verify_success(row):
    model, receipt = row['model'], row['receipt']
    if set(row) != {'model', 'status', 'receipt_path', 'receipt_sha256', 'receipt'} or row['status'] != 'retrieval_verified':
        raise ValueError('Invalid successful disposition schema')
    if receipt.get('status') != 'verified':
        raise ValueError('Unexpected receipt status')
    for field in ('model_id', 'version', 'sequence_sha256', 'length'):
        if receipt.get(field) != model.get(field):
            raise ValueError('Receipt identity differs from queued model')
    if receipt.get('url') != model.get('pae_url'):
        raise ValueError('Receipt URL differs from queued model')
    receipt_path = Path(row['receipt_path'])
    if digest(receipt_path) != row['receipt_sha256'] or json.loads(receipt_path.read_text()) != receipt:
        raise ValueError('Receipt file differs from disposition')
    matrix = Path(receipt['path'])
    compressed = matrix.read_bytes()
    if len(compressed) != receipt['compressed_bytes'] or digest(matrix) != receipt['gzip_sha256']:
        raise ValueError('Compressed PAE differs from receipt')
    raw = gzip.decompress(compressed)
    if len(raw) != receipt['json_bytes'] or hashlib.sha256(raw).hexdigest() != receipt['json_sha256']:
        raise ValueError('PAE JSON differs from receipt')
    return receipt['compressed_bytes'], receipt['json_bytes'], validate_matrix(raw, model)


def rows(path, tolerate_unfinalized=False):
    try:
        with gzip.open(path, 'rt', encoding='utf-8') as handle:
            for line in handle:
                if not line.endswith('\n'):
                    raise ValueError('Partial JSON row in disposition archive')
                yield json.loads(line)
    except EOFError:
        if not tolerate_unfinalized:
            raise


def audit(queue, old_manifest, recovery_manifest, old_preparation, recovery_receipt, output, receipt, script):
    if output.exists() or receipt.exists():
        raise FileExistsError('Fresh output and receipt paths required')
    old = json.loads(old_preparation.read_text())
    recovery = json.loads(recovery_receipt.read_text())
    if old.get('status') != 'prepared_interrupted_full_atlas_pae_recovery_queue':
        raise ValueError('Recovery preparation receipt required')
    if recovery.get('status') != 'completed_full_atlas_missing_pae_recovery_attempts_pending_independent_readback':
        raise ValueError('Completed recovery producer receipt required')
    if digest(old_manifest) != old['retained_manifest_sha256']:
        raise ValueError('Preserved manifest digest differs from preparation receipt')
    output.mkdir(parents=True)
    database = output / 'source-key-index.sqlite3'
    connection = sqlite3.connect(database)
    try:
        connection.execute('PRAGMA journal_mode=OFF')
        connection.execute('PRAGMA synchronous=OFF')
        connection.execute('CREATE TABLE source (k TEXT PRIMARY KEY, model TEXT NOT NULL, state TEXT NOT NULL)')
        expected = 0
        with queue.open() as handle:
            for raw in handle:
                model = json.loads(raw)
                connection.execute('INSERT INTO source VALUES (?, ?, ?)', (model_key(model), json.dumps(model, separators=(',', ':')), 'missing'))
                expected += 1
                if expected % 10000 == 0:
                    connection.commit()
        connection.commit()
        totals, errors = Counter(), Counter()
        compressed = json_bytes = entries = 0

        def check_manifest(path, label, unfinalized):
            nonlocal compressed, json_bytes, entries
            for row in rows(path, tolerate_unfinalized=unfinalized):
                model, status = row.get('model'), row.get('status')
                if not isinstance(model, dict) or status not in {'retrieval_verified', 'retrieval_failed'}:
                    raise ValueError('Invalid disposition')
                found = connection.execute('SELECT model, state FROM source WHERE k = ?', (model_key(model),)).fetchone()
                if found is None or json.loads(found[0]) != model:
                    raise ValueError('Disposition model differs from source queue')
                if found[1] != 'missing':
                    raise ValueError('Duplicate model disposition across manifests')
                if status == 'retrieval_verified':
                    cbytes, rbytes, cells = verify_success(row)
                    compressed += cbytes; json_bytes += rbytes; entries += cells
                else:
                    if set(row) != {'model', 'status', 'error_type', 'error'} or not isinstance(row['error_type'], str):
                        raise ValueError('Invalid failure disposition schema')
                    errors[row['error_type']] += 1
                connection.execute('UPDATE source SET state = ? WHERE k = ?', (label, model_key(model)))
                totals[label + '_' + status] += 1
                if sum(totals.values()) % 1000 == 0:
                    connection.commit()
        check_manifest(old_manifest, 'preserved', True)
        check_manifest(recovery_manifest, 'recovery', False)
        connection.commit()
        missing = connection.execute("SELECT COUNT(*) FROM source WHERE state = 'missing'").fetchone()[0]
        if missing:
            raise ValueError(f'{missing} source models lack a disposition')
        if totals['preserved_retrieval_verified'] != old['retained_dispositions']['retrieval_verified'] or totals['preserved_retrieval_failed'] != old['retained_dispositions']['retrieval_failed']:
            raise ValueError('Preserved disposition census differs from preparation receipt')
        if totals['recovery_retrieval_verified'] != recovery['counts']['verified'] or totals['recovery_retrieval_failed'] != recovery['counts']['failed']:
            raise ValueError('Recovery disposition census differs from producer receipt')
        if sum(totals.values()) != expected:
            raise ValueError('Combined disposition census differs from source queue')
        result = dict(status='passed_order_independent_full_atlas_pae_recovery_readback',
                      checked_utc=datetime.now(timezone.utc).isoformat(), expected_models=expected,
                      dispositions=dict(totals), error_types=dict(errors),
                      verified_compressed_bytes=compressed, verified_json_bytes=json_bytes,
                      verified_directional_entries=entries, source_queue_sha256=digest(queue),
                      preserved_manifest_sha256=digest(old_manifest), recovery_manifest_sha256=digest(recovery_manifest),
                      recovery_receipt_sha256=digest(recovery_receipt), independent_reader_sha256=digest(script),
                      scientific_eligibility=False, gpu=False, new_predictions=0,
                      scope='Independent source-key audit of every preserved and recovery disposition. Successful rows require exact source identity, retained receipt and gzip/JSON hashes, square finite nonnegative matrices and the AFDB export bound. Failed rows remain explicit. This is not confidence calibration, structural clustering, homology or evolutionary inference.')
        with receipt.open('x') as handle:
            json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
        (output / 'state.json').write_text(json.dumps(dict(stage='complete', dispositions=dict(totals)), indent=2) + '\n')
        return result
    finally:
        connection.close()


def self_test(script):
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary); queue = root / 'queue.jsonl'; old = root / 'old.jsonl.gz'; recovery = root / 'recovery.jsonl.gz'
        models = [dict(model_id=f'M{i}', version=6, sequence_sha256=str(i), length=2, pae_url=f'https://x/{i}') for i in range(3)]
        queue.write_text(''.join(json.dumps(model) + '\n' for model in models))
        with gzip.open(old, 'wt') as handle:
            handle.write(json.dumps(dict(model=models[1], status='retrieval_failed', error_type='HTTPError', error='429')) + '\n')
        with gzip.open(recovery, 'wt') as handle:
            for model in (models[2], models[0]):
                handle.write(json.dumps(dict(model=model, status='retrieval_failed', error_type='HTTPError', error='429')) + '\n')
        preparation = root / 'preparation.json'; preparation.write_text(json.dumps(dict(status='prepared_interrupted_full_atlas_pae_recovery_queue', retained_manifest_sha256=digest(old), retained_dispositions=dict(retrieval_verified=0, retrieval_failed=1))))
        producer = root / 'producer.json'; producer.write_text(json.dumps(dict(status='completed_full_atlas_missing_pae_recovery_attempts_pending_independent_readback', counts=dict(verified=0, failed=2))))
        result = audit(queue, old, recovery, preparation, producer, root / 'out', root / 'receipt.json', script)
        assert result['expected_models'] == 3 and result['dispositions']['preserved_retrieval_failed'] == 1 and result['dispositions']['recovery_retrieval_failed'] == 2
    return 'out_of_order_cross_manifest_source_key_audit'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--queue', type=Path); parser.add_argument('--preserved-manifest', type=Path)
    parser.add_argument('--recovery-manifest', type=Path); parser.add_argument('--preparation', type=Path)
    parser.add_argument('--recovery-receipt', type=Path); parser.add_argument('--output', type=Path); parser.add_argument('--receipt', type=Path)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(dict(status='passed_order_independent_pae_reader_literal_control', check=self_test(Path(__file__))), indent=2)); return
    required = (args.queue, args.preserved_manifest, args.recovery_manifest, args.preparation, args.recovery_receipt, args.output, args.receipt)
    if not all(required):
        parser.error('all source, manifest, receipt and output paths are required')
    print(json.dumps(audit(args.queue, args.preserved_manifest, args.recovery_manifest, args.preparation, args.recovery_receipt, args.output, args.receipt, Path(__file__)), indent=2))


if __name__ == '__main__':
    main()
