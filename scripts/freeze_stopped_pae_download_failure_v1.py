#!/usr/bin/env python3
"""Preserve a stopped V1 retrieval prefix and unreceipted bytes without adopting them."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import zlib

import psutil

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); assert not args.receipt.exists()
    review_path = Path('metadata/full_atlas_missing_pae_stop_review_20261005_v1.json')
    review = json.loads(review_path.read_text()); launch = review['launch']
    assert review['stop_returncode'] == 0 and review['original_wrapper_still_live'] is False
    try:
        p = psutil.Process(launch['pid'])
        assert p.create_time() != launch['created'] or p.cmdline() != launch['cmdline']
    except psutil.NoSuchProcess:
        pass
    plan_path = Path(launch['plan']); assert sha(plan_path) == launch['plan_sha256']
    plan = json.loads(plan_path.read_text()); pins = dict(plan['pins']); verify(pins)
    source = Path(plan['output'])/'retrieval_dispositions.jsonl.gz'
    blob = source.read_bytes(); decoder = zlib.decompressobj(wbits=31); prefix = decoder.decompress(blob)
    newline = prefix.rfind(b'\n'); complete = prefix[:newline+1]; trailing = prefix[newline+1:]
    rows = [json.loads(line) for line in complete.splitlines()]
    assert rows and not decoder.unused_data and len(rows) >= review['before_state']['counts']['attempted']
    models = {}
    with Path(plan['queue']).open() as handle:
        for index, line in enumerate(handle):
            if index >= len(rows)+plan['future_window']: break
            model = json.loads(line); models[model['model_id']+'-v'+str(model['version'])] = model
    seen = set(); errors = Counter()
    for row in rows:
        model = row['model']; key = model['model_id']+'-v'+str(model['version'])
        assert key not in seen and models[key] == model; seen.add(key)
        assert row['status'] == 'retrieval_failed' and row['error_type'] == 'ValueError'
        assert 'one path is relative and the other is absolute' in row['error']
        errors[row['error_type']] += 1
    root = Path('results/structures/full-atlas-stopped-pae-failure-review-20261005-v1')
    assert not root.exists(); root.mkdir(parents=True)
    recovered = root/'complete_original_error_prefix.jsonl'
    recovered.write_bytes(complete)
    orphan_manifest = root/'unreceipted_byte_inventory.jsonl'
    count = 0; size = 0; partial_count = 0
    cache = Path(plan['cache']); assert not list(cache.rglob('*.receipt.json'))
    with orphan_manifest.open('x') as handle:
        for path in sorted(cache.rglob('*')):
            if not path.is_file(): continue
            match = re.fullmatch(r'(.+-v[0-9]+)\.json\.(gz|tmp)', path.name)
            assert match and match[1] in models
            row = dict(path=str(path), sha256=sha(path), bytes=path.stat().st_size,
                       model=models[match[1]], status='unreceipted_bytes_not_admitted',
                       partial_file=match[2]=='tmp', body_content_validated=False)
            handle.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
            count += 1; size += row['bytes']; partial_count += int(row['partial_file'])
    assert count >= len(rows) and count <= len(rows)+plan['future_window']
    for path in (source, recovered, orphan_manifest, review_path, plan_path, Path(__file__), Path(review['journal_path'])):
        bind(pins, path)
    verify(pins)
    result = dict(status='preserved_stopped_original_pae_error_prefix_and_unreceipted_inventory',
                  checked_utc=datetime.now(timezone.utc).isoformat(), original_invocation=launch['invocation_id'],
                  original_tool_session_id=launch['original_tool_session_id'], original_tool_terminal_exit_code=None,
                  native_success_claimed=False, recovered_complete_error_rows=len(rows), error_types=dict(errors),
                  gzip_end_of_stream_present=decoder.eof, trailing_decompressed_bytes=len(trailing),
                  trailing_decompressed_sha256=hashlib.sha256(trailing).hexdigest(),
                  unreceipted_files=count, unreceipted_bytes=size, partial_files=partial_count,
                  accepted_cached_matrices=0, source_hashes=pins, scientific_eligibility=False,
                  scope='Stopped exact original V1; original gzip remains unchanged, including any missing trailer. '
                        'Only complete original JSON lines are recovered into a new artifact and matched to '
                        'original bounded queue prefix; all unreceipted files are byte-inventoried. No file '
                        'is accepted as a matrix, no HTTP provenance or native zero exit is invented, no '
                        'original command is repeated. V2 retrieves full original queue in a separate cache.')
    with args.receipt.open('x') as handle: json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}, indent=2))


if __name__ == '__main__': main()
