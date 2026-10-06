#!/usr/bin/env python3
"""Read-only integrity checkpoint for a completed prefix of a live PAE queue."""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def validate(raw, model):
    payload = json.loads(raw)
    assert isinstance(payload, list) and len(payload) == 1 and isinstance(payload[0], dict)
    matrix = payload[0]['predicted_aligned_error']; maximum = payload[0]['max_predicted_aligned_error']
    assert isinstance(matrix, list) and len(matrix) == model['length']
    assert isinstance(maximum, (int, float)) and math.isfinite(float(maximum)) and maximum > 0
    entries = 0
    for row in matrix:
        assert isinstance(row, list) and len(row) == model['length']
        for value in row:
            assert isinstance(value, (int, float)) and math.isfinite(float(value)) and 0 <= value <= maximum + .51
            entries += 1
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--queue', type=Path, required=True)
    parser.add_argument('--dispositions', type=Path, required=True)
    parser.add_argument('--records', type=int, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert args.records > 0 and not args.receipt.exists()
    compressed = json_bytes = entries = 0
    with gzip.GzipFile(args.dispositions, 'rb') as stream:
        rows = [json.loads(stream.readline()) for _ in range(args.records)]
    keys = {(row['model']['model_id'], row['model']['version']) for row in rows}
    assert len(keys) == args.records
    queued_models = {}
    with args.queue.open() as queued:
        for line in queued:
            model = json.loads(line); key = (model['model_id'], model['version'])
            if key in keys:
                queued_models[key] = model
                if len(queued_models) == len(keys):
                    break
    assert set(queued_models) == keys
    for row in rows:
            model = queued_models[(row['model']['model_id'], row['model']['version'])]
            assert row['status'] == 'retrieval_verified' and row['model'] == model
            receipt = row['receipt']; receipt_path = Path(row['receipt_path']); matrix_path = Path(receipt['path'])
            assert sha(receipt_path) == row['receipt_sha256'] and json.loads(receipt_path.read_text()) == receipt
            for key in ('model_id', 'version', 'sequence_sha256', 'length'):
                assert receipt[key] == model[key]
            assert receipt['url'] == model['pae_url']
            body = matrix_path.read_bytes(); raw = gzip.decompress(body)
            assert len(body) == receipt['compressed_bytes'] and sha(matrix_path) == receipt['gzip_sha256']
            assert len(raw) == receipt['json_bytes'] and hashlib.sha256(raw).hexdigest() == receipt['json_sha256']
            compressed += len(body); json_bytes += len(raw); entries += validate(raw, model)
    result = dict(status='passed_live_missing_pae_prefix_integrity_checkpoint', records=args.records,
                  queue_sha256=sha(args.queue), disposition_stream_bytes=args.dispositions.stat().st_size,
                  verified_compressed_bytes=compressed, verified_json_bytes=json_bytes,
                  verified_directional_entries=entries, script_sha256=sha(Path(__file__)),
                  scope='Read-only validation of a completed prefix while the gzip disposition stream remains live. Each completion-order row is matched back to its unique original queue record, then its receipt, compressed/decompressed bytes, matrix dimensions and finite nonnegative values are checked. This checkpoint does not read the trailing stream end, establish whole-queue completion, alter cache files, or make confidence/biological inferences.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
