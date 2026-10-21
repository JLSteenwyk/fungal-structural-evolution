#!/usr/bin/env python3
"""Make an ordered recovery queue from an interrupted full-atlas PAE retrieval.

The downloader records rows as HTTP futures finish, which is intentionally not
the source-queue order.  This helper preserves that fact: it validates every
complete retained row against the immutable source queue and emits only models
without a retained verified disposition (both prior failures and unattempted
models).  It does no network I/O and never changes the original cache.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def key(model):
    return f"{model['model_id']}-v{model['version']}"


def read_retained_rows(path):
    """Yield every newline-complete gzip row, tolerating a missing final trailer."""
    try:
        with gzip.open(path, 'rt', encoding='utf-8') as handle:
            for line in handle:
                if not line.endswith('\n'):
                    raise ValueError('Interrupted manifest ends with a partial JSON row')
                yield json.loads(line)
    except EOFError:
        # The original writer was interrupted before gzip.finalize().  Every
        # line already yielded was complete and is checked below; no trailing
        # bytes are admitted by this recovery helper.
        return


def prepare(queue, retained, output, receipt, script):
    if output.exists() or receipt.exists():
        raise FileExistsError('Fresh output directory and receipt path required')
    output.mkdir(parents=True)
    database = output / 'source-index.sqlite3'
    connection = sqlite3.connect(database)
    try:
        connection.execute('PRAGMA journal_mode=OFF')
        connection.execute('PRAGMA synchronous=OFF')
        connection.execute('CREATE TABLE source (k TEXT PRIMARY KEY, ordinal INTEGER UNIQUE, raw TEXT, state TEXT NOT NULL)')
        expected = 0
        with queue.open(encoding='utf-8') as handle:
            for ordinal, raw in enumerate(handle):
                model = json.loads(raw)
                connection.execute('INSERT INTO source VALUES (?, ?, ?, ?)', (key(model), ordinal, raw, 'unattempted'))
                expected += 1
                if ordinal % 10000 == 0:
                    connection.commit()
        connection.commit()
        observed = Counter()
        for row in read_retained_rows(retained):
            if set(row) - {'model', 'status', 'receipt_path', 'receipt_sha256', 'receipt', 'error_type', 'error'}:
                raise ValueError('Unexpected retained disposition fields')
            model, status = row.get('model'), row.get('status')
            if not isinstance(model, dict) or status not in {'retrieval_verified', 'retrieval_failed'}:
                raise ValueError('Invalid retained disposition')
            matched = connection.execute('SELECT raw, state FROM source WHERE k = ?', (key(model),)).fetchone()
            if matched is None or json.loads(matched[0]) != model:
                raise ValueError('Retained model is absent from or differs from source queue')
            if matched[1] != 'unattempted':
                raise ValueError('Duplicate retained model disposition')
            if status == 'retrieval_verified':
                required = {'model', 'status', 'receipt_path', 'receipt_sha256', 'receipt'}
            else:
                required = {'model', 'status', 'error_type', 'error'}
            if set(row) != required:
                raise ValueError('Disposition schema differs from its status')
            connection.execute('UPDATE source SET state = ? WHERE k = ?',
                               ('verified' if status == 'retrieval_verified' else 'failed', key(model)))
            observed[status] += 1
            if sum(observed.values()) % 10000 == 0:
                connection.commit()
        connection.commit()
        states = dict(connection.execute('SELECT state, COUNT(*) FROM source GROUP BY state'))
        if states.get('verified', 0) != observed['retrieval_verified'] or states.get('failed', 0) != observed['retrieval_failed']:
            raise ValueError('Retained disposition census differs from indexed states')
        recovery_queue = output / 'recovery_queue.jsonl'
        with recovery_queue.open('x', encoding='utf-8') as handle:
            for (raw,) in connection.execute("SELECT raw FROM source WHERE state != 'verified' ORDER BY ordinal"):
                handle.write(raw if raw.endswith('\n') else raw + '\n')
        recovery_models = expected - states.get('verified', 0)
        if sum(1 for _ in recovery_queue.open(encoding='utf-8')) != recovery_models:
            raise ValueError('Recovery queue count differs from source-index census')
        result = {
            'status': 'prepared_interrupted_full_atlas_pae_recovery_queue',
            'checked_utc': datetime.now(timezone.utc).isoformat(),
            'source_queue': str(queue), 'source_queue_sha256': sha256(queue),
            'retained_manifest': str(retained), 'retained_manifest_sha256': sha256(retained),
            'expected_models': expected, 'retained_dispositions': dict(observed),
            'source_states': states, 'recovery_models': recovery_models,
            'recovery_queue': str(recovery_queue), 'recovery_queue_sha256': sha256(recovery_queue),
            'script_sha256': sha256(script),
            'scope': 'Ordered recovery inventory only. Complete retained rows are source-bound and duplicate-checked; '
                     'the queue contains prior failures and unattempted source models. It neither validates retained '
                     'PAE matrices nor retrieves, deletes, accepts or modifies any PAE file.'
        }
        with receipt.open('x', encoding='utf-8') as handle:
            json.dump(result, handle, indent=2, allow_nan=False)
            handle.write('\n')
        return result
    finally:
        connection.close()


def self_test(script):
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        queue = root / 'queue.jsonl'
        models = [dict(model_id=f'M{i}', version=6, sequence_sha256=str(i), length=i + 2, pae_url=f'https://x/{i}') for i in range(3)]
        queue.write_text(''.join(json.dumps(row) + '\n' for row in models))
        retained = root / 'rows.jsonl.gz'
        with gzip.open(retained, 'wt', encoding='utf-8') as handle:
            handle.write(json.dumps(dict(model=models[1], status='retrieval_verified', receipt_path='r', receipt_sha256='s', receipt={})) + '\n')
            handle.write(json.dumps(dict(model=models[0], status='retrieval_failed', error_type='HTTPError', error='429')) + '\n')
        output, receipt = root / 'output', root / 'receipt.json'
        result = prepare(queue, retained, output, receipt, script)
        assert result['recovery_models'] == 2
        assert [json.loads(row)['model_id'] for row in (output / 'recovery_queue.jsonl').read_text().splitlines()] == ['M0', 'M2']
    return 'source-order-independent_recovery_queue'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--queue', type=Path)
    parser.add_argument('--retained-manifest', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--receipt', type=Path)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps({'status': 'passed_recovery_queue_literal_control', 'check': self_test(Path(__file__))}, indent=2))
        return
    if not all((args.queue, args.retained_manifest, args.output, args.receipt)):
        parser.error('queue, retained manifest, output, and receipt are required')
    print(json.dumps(prepare(args.queue, args.retained_manifest, args.output, args.receipt, Path(__file__)), indent=2))


if __name__ == '__main__':
    main()
