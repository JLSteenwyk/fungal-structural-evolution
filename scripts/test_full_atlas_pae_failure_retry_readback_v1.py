#!/usr/bin/env python3
"""Exercise successful-matrix validation in the immutable PAE retry closure reader."""
import gzip
import hashlib
import json
from pathlib import Path
import tempfile

from readback_full_atlas_pae_failure_retry_v1 import audit, digest


def main():
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        queue, old, recovery, retried = (root / name for name in ('queue.jsonl', 'old.jsonl.gz', 'recovery.jsonl.gz', 'retry.jsonl.gz'))
        model = dict(model_id='M1', version=6, sequence_sha256='abc', length=2, pae_url='https://x/M1')
        queue.write_text(json.dumps(model) + '\n')
        failed = dict(model=model, status='retrieval_failed', error_type='HTTPError', error='429')
        for path, records in ((old, [failed]), (recovery, [])):
            with gzip.open(path, 'wt') as handle:
                for row in records:
                    handle.write(json.dumps(row) + '\n')
        raw = json.dumps([dict(predicted_aligned_error=[[0.0, 1.0], [1.0, 0.0]], max_predicted_aligned_error=1.0)], separators=(',', ':')).encode()
        matrix = root / 'matrix.json.gz'; matrix.write_bytes(gzip.compress(raw, mtime=0))
        matrix_receipt = dict(model_id='M1', version=6, sequence_sha256='abc', length=2, status='verified', url=model['pae_url'], path=str(matrix), gzip_sha256=digest(matrix), json_sha256=hashlib.sha256(raw).hexdigest(), compressed_bytes=matrix.stat().st_size, json_bytes=len(raw))
        receipt_path = root / 'matrix.receipt.json'; receipt_path.write_text(json.dumps(matrix_receipt))
        success = dict(model=model, status='retrieval_verified', receipt_path=str(receipt_path), receipt_sha256=digest(receipt_path), receipt=matrix_receipt)
        with gzip.open(retried, 'wt') as handle:
            handle.write(json.dumps(success) + '\n')
        base = root / 'base.json'; base.write_text(json.dumps(dict(status='passed_order_independent_full_atlas_pae_recovery_readback', dispositions=dict(preserved_retrieval_failed=1, recovery_retrieval_failed=0))))
        producer = root / 'producer.json'; producer.write_text(json.dumps(dict(status='completed_full_atlas_pae_explicit_failure_retry_pending_independent_readback', retry_manifest_sha256=digest(retried), expected_models=1, counts=dict(attempted=1, verified=1, failed=0))))
        result = audit(queue, old, recovery, retried, base, producer, root / 'out', root / 'closure.json', Path(__file__).with_name('readback_full_atlas_pae_failure_retry_v1.py'))
        assert result['retained_base_failures'] == 1
        assert result['final_explicit_failures'] == 0
        assert result['retry_verified_directional_entries'] == 4
    print(json.dumps(dict(status='passed_successful_pae_retry_closure_fixture', matrix_entries=4), indent=2))


if __name__ == '__main__':
    main()
