#!/usr/bin/env python3
"""Build an ordered retry queue from complete PAE disposition manifests."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def key(model):
    return f"{model['model_id']}-v{model['version']}"


def rows(path):
    with gzip.open(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if not line.endswith('\n'):
                raise ValueError('Partial disposition row')
            yield json.loads(line)


def prepare(queue, manifests, output, receipt, script):
    if output.exists() or receipt.exists():
        raise FileExistsError('Fresh output and receipt required')
    output.mkdir(parents=True)
    index = output / 'source-index.sqlite3'
    db = sqlite3.connect(index)
    try:
        db.execute('PRAGMA journal_mode=OFF'); db.execute('PRAGMA synchronous=OFF')
        db.execute('CREATE TABLE source (k TEXT PRIMARY KEY, ordinal INTEGER UNIQUE, raw TEXT, state TEXT NOT NULL)')
        expected = 0
        with queue.open() as handle:
            for ordinal, raw in enumerate(handle):
                model = json.loads(raw)
                db.execute('INSERT INTO source VALUES (?, ?, ?, ?)', (key(model), ordinal, raw, 'unseen'))
                expected += 1
        db.commit(); counts = Counter(); errors = Counter()
        for manifest in manifests:
            for row in rows(manifest):
                model, status = row.get('model'), row.get('status')
                if not isinstance(model, dict) or status not in {'retrieval_verified', 'retrieval_failed'}:
                    raise ValueError('Invalid disposition')
                found = db.execute('SELECT raw, state FROM source WHERE k=?', (key(model),)).fetchone()
                if found is None or json.loads(found[0]) != model:
                    raise ValueError('Disposition differs from source queue')
                if found[1] != 'unseen':
                    raise ValueError('Duplicate disposition across manifests')
                if status == 'retrieval_failed':
                    if set(row) != {'model', 'status', 'error_type', 'error'}:
                        raise ValueError('Invalid failure schema')
                    state = 'failed'; errors[row['error_type']] += 1
                else:
                    if set(row) != {'model', 'status', 'receipt_path', 'receipt_sha256', 'receipt'}:
                        raise ValueError('Invalid success schema')
                    state = 'verified'
                db.execute('UPDATE source SET state=? WHERE k=?', (state, key(model))); counts[state] += 1
        db.commit()
        unseen = db.execute("SELECT COUNT(*) FROM source WHERE state='unseen'").fetchone()[0]
        if unseen:
            raise ValueError(f'{unseen} source models lack a disposition')
        retry = output / 'retry_queue.jsonl'
        with retry.open('x') as handle:
            for raw, in db.execute("SELECT raw FROM source WHERE state='failed' ORDER BY ordinal"):
                handle.write(raw if raw.endswith('\n') else raw + '\n')
        result = dict(status='prepared_full_atlas_pae_explicit_failure_retry_queue', checked_utc=datetime.now(timezone.utc).isoformat(),
                      expected_models=expected, dispositions=dict(counts), error_types=dict(errors), retry_models=counts['failed'],
                      source_queue=str(queue), source_queue_sha256=digest(queue), manifests=[dict(path=str(p), sha256=digest(p)) for p in manifests],
                      retry_queue=str(retry), retry_queue_sha256=digest(retry), script_sha256=digest(script),
                      scope='Read-only source-key reconciliation of completed disposition manifests. The queue contains only explicit retained failures; it neither retrieves nor validates matrices.')
        with receipt.open('x') as handle: json.dump(result, handle, indent=2); handle.write('\n')
        return result
    finally:
        db.close()


def self_test(script):
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp); queue = root/'queue'; first = root/'first.gz'; second = root/'second.gz'
        models = [dict(model_id=f'M{i}', version=6) for i in range(3)]
        queue.write_text(''.join(json.dumps(x)+'\n' for x in models))
        with gzip.open(first, 'wt') as h:
            h.write(json.dumps(dict(model=models[1], status='retrieval_failed', error_type='HTTPError', error='429'))+'\n')
            h.write(json.dumps(dict(model=models[0], status='retrieval_verified', receipt_path='x', receipt_sha256='x', receipt={}))+'\n')
        with gzip.open(second, 'wt') as h:
            h.write(json.dumps(dict(model=models[2], status='retrieval_verified', receipt_path='y', receipt_sha256='y', receipt={}))+'\n')
        r = prepare(queue, [first, second], root/'out', root/'receipt', script)
        assert r['retry_models'] == 1
        assert json.loads((root/'out/retry_queue.jsonl').read_text())['model_id'] == 'M1'
    return 'cross_manifest_explicit_failure_only_queue'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--queue', type=Path); p.add_argument('--manifest', type=Path, action='append')
    p.add_argument('--output', type=Path); p.add_argument('--receipt', type=Path); p.add_argument('--self-test', action='store_true')
    a = p.parse_args()
    if a.self_test:
        print(json.dumps(dict(status='passed_failure_retry_queue_literal_control', check=self_test(Path(__file__)), indent=2))); return
    if not all((a.queue, a.manifest, a.output, a.receipt)):
        p.error('queue, one or more manifests, output and receipt are required')
    print(json.dumps(prepare(a.queue, a.manifest, a.output, a.receipt, Path(__file__)), indent=2))


if __name__ == '__main__':
    main()
