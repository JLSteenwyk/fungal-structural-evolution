#!/usr/bin/env python3
"""Independently replay a completed full missing-AFDB-PAE retrieval queue."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import tempfile


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def validate_matrix(raw, model):
    """Validate the published AFDB monomer PAE JSON without retriever imports."""
    payload = json.loads(raw)
    if not isinstance(payload, list) or len(payload) != 1 or not isinstance(payload[0], dict):
        raise ValueError('Invalid PAE JSON container')
    body = payload[0]
    values, maximum = body.get('predicted_aligned_error'), body.get('max_predicted_aligned_error')
    length = model['length']
    if not isinstance(values, list) or len(values) != length or not isinstance(maximum, (int, float)):
        raise ValueError('Invalid PAE dimensions or declared maximum')
    if not math.isfinite(float(maximum)) or float(maximum) <= 0:
        raise ValueError('Invalid declared PAE maximum')
    count = 0
    for row in values:
        if not isinstance(row, list) or len(row) != length:
            raise ValueError('Non-square PAE matrix')
        for value in row:
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)) or value < 0:
                raise ValueError('Non-finite or negative PAE value')
            # Preserve the source's established AFDB export-rounding allowance.
            if value > maximum + .51:
                raise ValueError('PAE value exceeds declared export bound')
            count += 1
    return count


def verify_success(row):
    model, receipt = row['model'], row['receipt']
    if row['status'] != 'retrieval_verified' or receipt.get('status') != 'verified':
        raise ValueError('Unexpected success disposition')
    for field in ('model_id', 'version', 'sequence_sha256', 'length'):
        if receipt.get(field) != model.get(field):
            raise ValueError('Receipt identity differs from queued model')
    if receipt.get('url') != model.get('pae_url'):
        raise ValueError('Receipt URL differs from queued model')
    receipt_path = Path(row['receipt_path'])
    if digest(receipt_path) != row['receipt_sha256'] or json.loads(receipt_path.read_text()) != receipt:
        raise ValueError('Receipt file differs from disposition')
    matrix_path = Path(receipt['path'])
    compressed = matrix_path.read_bytes()
    if len(compressed) != receipt['compressed_bytes'] or digest(matrix_path) != receipt['gzip_sha256']:
        raise ValueError('Compressed PAE file differs from receipt')
    raw = gzip.decompress(compressed)
    if len(raw) != receipt['json_bytes'] or hashlib.sha256(raw).hexdigest() != receipt['json_sha256']:
        raise ValueError('PAE JSON differs from receipt')
    entries = validate_matrix(raw, model)
    return receipt['compressed_bytes'], receipt['json_bytes'], entries


def self_test():
    model = dict(model_id='OFFLINE', version=6, sequence_sha256='x', length=2,
                 pae_url='https://example.invalid/pae.json')
    payload = json.dumps([dict(predicted_aligned_error=[[0, 3.5], [2, 0]],
                              max_predicted_aligned_error=4)]).encode()
    assert validate_matrix(payload, model) == 4
    checks = ['literal_square_finite_nonnegative_matrix']
    for altered in (
        b'[]',
        json.dumps([dict(predicted_aligned_error=[[0], [2]], max_predicted_aligned_error=4)]).encode(),
        json.dumps([dict(predicted_aligned_error=[[0, -1], [2, 0]], max_predicted_aligned_error=4)]).encode(),
        json.dumps([dict(predicted_aligned_error=[[0, 5], [2, 0]], max_predicted_aligned_error=4)]).encode(),
    ):
        try:
            validate_matrix(altered, model)
        except ValueError:
            checks.append('reject_malformed_or_invalid_value')
        else:
            raise AssertionError('Invalid PAE fixture accepted')
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary); matrix = root/'matrix.json.gz'; receipt_path = root/'receipt.json'
        with matrix.open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as handle:
            handle.write(payload)
        # gzip.open does not expose mtime on every supported Python version; the
        # payload, not its transport timestamp, is what the receipt binds.
        receipt = dict(**model, status='verified', url=model['pae_url'], path=str(matrix),
                       gzip_sha256=digest(matrix), json_sha256=hashlib.sha256(payload).hexdigest(),
                       compressed_bytes=matrix.stat().st_size, json_bytes=len(payload))
        receipt_path.write_text(json.dumps(receipt))
        row = dict(model=model, status='retrieval_verified', receipt_path=str(receipt_path),
                   receipt_sha256=digest(receipt_path), receipt=receipt)
        assert verify_success(row)[2] == 4
        checks.append('literal_receipt_and_gzip_hash_binding')
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer-receipt', type=Path)
    parser.add_argument('--producer-plan', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--receipt', type=Path)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    checks = self_test()
    if args.self_test:
        print(json.dumps(dict(status='passed_full_missing_pae_readback_literal_controls', checks=checks), indent=2))
        return
    if not all((args.producer_receipt, args.producer_plan, args.output, args.receipt)):
        parser.error('producer paths, output and receipt are required unless --self-test is used')
    if args.output.exists() or args.receipt.exists():
        raise FileExistsError('Fresh output and receipt paths required')
    plan, producer = (json.loads(path.read_text()) for path in (args.producer_plan, args.producer_receipt))
    if producer.get('status') != 'completed_complete_missing_afdb_pae_retrieval_attempts_pending_independent_readback':
        raise ValueError('Completed retrieval receipt required')
    if producer.get('source_hashes', {}).get(str(args.producer_plan)) != digest(args.producer_plan):
        raise ValueError('Producer plan digest differs')
    expected = plan['expected_models']; queue = Path(plan['queue']); dispositions = Path(plan['output'])/'retrieval_dispositions.jsonl.gz'
    counts, errors = Counter(), Counter(); compressed = json_bytes = entries = 0
    args.output.mkdir(parents=True)
    with queue.open() as queued, gzip.open(dispositions, 'rt') as observed:
        for number, line in enumerate(queued, 1):
            model, row = json.loads(line), json.loads(next(observed))
            if row.get('model') != model:
                raise ValueError('Disposition order or source record differs')
            counts[row.get('status')] += 1
            if row.get('status') == 'retrieval_verified':
                cbytes, rbytes, cells = verify_success(row)
                compressed += cbytes; json_bytes += rbytes; entries += cells
            elif row.get('status') == 'retrieval_failed':
                if set(row) != {'model', 'status', 'error_type', 'error'} or not isinstance(row['error_type'], str):
                    raise ValueError('Invalid retained retrieval-failure row')
                errors[row['error_type']] += 1
            else:
                raise ValueError('Unknown retrieval disposition')
            if number % 1000 == 0:
                (args.output/'state.json').write_text(json.dumps(dict(stage='independent_full_missing_pae_readback',
                    models=number, expected_models=expected, dispositions=dict(counts)), indent=2)+'\n')
        if next(observed, None) is not None:
            raise ValueError('Extra retrieval disposition rows')
    expected_dispositions = dict(retrieval_verified=producer['counts']['verified'],
                                 retrieval_failed=producer['counts']['failed'])
    observed_dispositions = dict(retrieval_verified=counts['retrieval_verified'],
                                 retrieval_failed=counts['retrieval_failed'])
    if sum(counts.values()) != expected or observed_dispositions != expected_dispositions:
        raise ValueError('Complete disposition census differs from producer')
    if dict(errors) != producer['error_types']:
        raise ValueError('Error-type census differs from producer')
    result = dict(status='passed_independent_full_missing_afdb_pae_retrieval_readback',
                  checked_utc=datetime.now(timezone.utc).isoformat(), checks=checks,
                  expected_models=expected, dispositions=dict(counts), error_types=dict(errors),
                  verified_compressed_bytes=compressed, verified_json_bytes=json_bytes,
                  verified_directional_entries=entries, producer_receipt_sha256=digest(args.producer_receipt),
                  producer_plan_sha256=digest(args.producer_plan), independent_reader_sha256=digest(Path(__file__)),
                  scientific_eligibility=False, gpu=False, new_predictions=0,
                  scope='Independent ordered replay of every completed missing-AFDB PAE queue row. Successful rows require exact source identity, receipt and gzip/JSON hashes, square finite nonnegative matrix values and the established AFDB export bound; failures remain explicit. No PAE symmetrization, confidence/context calibration, structure prediction, structural homology, or evolutionary inference.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    (args.output/'state.json').write_text(json.dumps(dict(stage='complete', dispositions=dict(counts)), indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
