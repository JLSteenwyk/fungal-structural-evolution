#!/usr/bin/env python3
"""Summarize retained PAE-recovery failures without altering the live archive.

The recovery journal is an append-only gzip stream. While its producer is
running, the last gzip member has no end marker yet. This reader reports all
complete JSON records available before that tail and calls the result a partial
observation instead of treating the expected EOF condition as data corruption.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any


DEFAULT_ARCHIVE = Path(
    'results/structures/full-atlas-missing-pae-recovery-20261021-v1/'
    'recovery_dispositions.jsonl.gz')
HTTP_STATUS = re.compile(r'HTTP Error (\d{3})\b')


def summarize(archive: Path) -> dict[str, Any]:
    """Read complete records in *archive*, tolerating a live gzip tail."""
    statuses: Counter[str] = Counter()
    error_types: Counter[str] = Counter()
    http_statuses: Counter[str] = Counter()
    records_read = 0
    incomplete_stream = False
    parse_error: str | None = None
    try:
        with gzip.open(archive, 'rt', encoding='utf-8') as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                try:
                    row = json.loads(line)
                except json.JSONDecodeError as exc:
                    parse_error = f'line {line_number}: {exc.msg}'
                    break
                if not isinstance(row, dict):
                    parse_error = f'line {line_number}: JSON record is not an object'
                    break
                records_read += 1
                status = row.get('status')
                if isinstance(status, str):
                    statuses[status] += 1
                error_type = row.get('error_type')
                if isinstance(error_type, str):
                    error_types[error_type] += 1
                error = row.get('error')
                if isinstance(error, str):
                    match = HTTP_STATUS.search(error)
                    if match:
                        http_statuses[match.group(1)] += 1
    except EOFError:
        incomplete_stream = True

    return {
        'status': 'pae_recovery_failure_summary_advisory',
        'archive': str(archive),
        'records_read_before_stream_end': records_read,
        'incomplete_live_gzip_tail': incomplete_stream,
        'json_parse_error': parse_error,
        'disposition_status_counts': dict(sorted(statuses.items())),
        'error_type_counts': dict(sorted(error_types.items())),
        'http_status_counts': dict(sorted(http_statuses.items())),
        'scope': ('Read-only operational observation. A live gzip tail can omit recent records; '
                  'this report neither retries models nor establishes PAE qualification or scientific eligibility.'),
    }


def self_test() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        archive = Path(temporary) / 'journal.jsonl.gz'
        rows = [
            {'status': 'verified'},
            {'status': 'retrieval_failed', 'error_type': 'HTTPError',
             'error': 'HTTP Error 429: Too Many Requests'},
        ]
        with gzip.open(archive, 'wt', encoding='utf-8') as handle:
            for row in rows:
                handle.write(json.dumps(row) + '\n')
        complete = summarize(archive)
        assert complete['records_read_before_stream_end'] == 2
        assert not complete['incomplete_live_gzip_tail']
        assert complete['http_status_counts'] == {'429': 1}
        archive.write_bytes(archive.read_bytes()[:-8])
        partial = summarize(archive)
        assert partial['records_read_before_stream_end'] == 2
        assert partial['incomplete_live_gzip_tail']


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    self_test()
    if args.self_test:
        print(json.dumps({'status': 'passed_pae_recovery_failure_summary_controls',
                          'checks': ['complete_gzip', 'live_gzip_tail', 'http_status_classification']}, indent=2))
        return
    print(json.dumps(summarize(args.archive), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
