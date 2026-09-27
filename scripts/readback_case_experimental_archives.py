#!/usr/bin/env python3
"""Independently stream-check every case-homolog archive after terminal retrieval."""
import gzip
import hashlib
import json
import subprocess
import time
from pathlib import Path

import psutil
from screen_duplication_domain_alignment_coverage import sha

ROOT = Path('results/experimental_structures/whole-domain-case-coordinates-20260927-v1')
LAUNCH = Path('metadata/case_experimental_coordinate_launch_20260927.json')
OUTPUT = Path('results/experimental_structures/whole-domain-case-coordinate-readback-20260927-v1')


def main():
    launch = json.loads(LAUNCH.read_text())
    sources = {str(LAUNCH): sha(LAUNCH)}
    while True:
        try:
            process = psutil.Process(launch['pid'])
            if process.create_time() != launch['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            if process.cmdline() != launch['cmdline']:
                raise ValueError('Producer command changed')
        except psutil.NoSuchProcess:
            break
        print('Waiting for exact coordinate producer', launch['pid'], flush=True)
        time.sleep(30)
    state = dict(line.split('=', 1) for line in subprocess.check_output(
        ['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState',
         '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    assert sha(launch['cmdline'][1]) == launch['script_sha256']
    config_path, receipt_path = ROOT / 'config.json', ROOT / 'receipt.json'
    config, receipt = (json.loads(p.read_text()) for p in (config_path, receipt_path))
    assert receipt['status'] == 'complete_frozen_experimental_coordinate_download'
    assert receipt['config_sha256'] == sha(config_path)
    metadata_root = Path('results/experimental_structures/whole-domain-case-metadata-20260927-v1')
    metadata_path = metadata_root / 'receipt.json'
    audit_path = Path('results/experimental_structures/whole-domain-case-subject-readback-20260927-v1/receipt.json')
    assert sha(metadata_path) == config['metadata_receipt_sha256']
    assert sha(audit_path) == config['subject_audit_receipt_sha256']
    for p in (config_path, receipt_path, metadata_path, audit_path):
        sources[str(p)] = sha(p)
    expected, excluded = set(), {}
    for item in json.loads(metadata_path.read_text())['responses']:
        if item['kind'] != 'entry':
            continue
        path = metadata_root / 'entry' / (item['identifier'] + '.json')
        assert sha(path) == item['response_sha256']
        sources[str(path)] = sha(path)
        methodology = json.loads(path.read_text())['rcsb_entry_info'].get('structure_determination_methodology', 'unknown')
        if methodology == 'experimental':
            expected.add(item['identifier'])
        else:
            excluded[item['identifier']] = methodology
    assert excluded == config['excluded_entries']
    assert sorted(expected) == config['entries']
    assert {r['entry_id'] for r in receipt['results']} == expected
    assert len(receipt['results']) == len(expected) == receipt['entries']
    checked = []
    for item in receipt['results']:
        entry = item['entry_id']
        path, rp = ROOT / (entry + '.cif.gz'), ROOT / (entry + '.receipt.json')
        r = json.loads(rp.read_text())
        assert sha(rp) == item['receipt_sha256']
        assert r['entry_id'] == entry and r['config_sha256'] == sha(config_path)
        assert r['url'] == 'https://files.rcsb.org/download/' + entry + '.cif.gz'
        assert sha(path) == r['gzip_sha256'] == item['gzip_sha256']
        assert path.stat().st_size == r['compressed_bytes'] == item['compressed_bytes']
        count, digest, first = 0, hashlib.sha256(), b''
        # Independent decompression to EOF verifies CRC and records plaintext hash.
        with gzip.open(path, 'rb') as handle:
            while chunk := handle.read(65536):
                if not first:
                    first = chunk
                digest.update(chunk)
                count += len(chunk)
        first_content = next(line.strip() for line in first.splitlines()
                             if line.strip() and not line.lstrip().startswith(b'#'))
        assert first_content.lower() == ('data_' + entry).lower().encode()
        assert count == r['uncompressed_bytes'] == item['uncompressed_bytes']
        sources[str(rp)] = sha(rp)
        checked.append(dict(entry_id=entry, gzip_sha256=item['gzip_sha256'],
                            uncompressed_sha256=digest.hexdigest(), uncompressed_bytes=count))
        print(len(checked), '/', len(expected), entry, flush=True)
    assert sum(x['uncompressed_bytes'] for x in checked) == receipt['uncompressed_bytes']
    assert sum(x['compressed_bytes'] for x in receipt['results']) == receipt['compressed_bytes']
    for path, digest in sources.items():
        assert sha(path) == digest
    OUTPUT.mkdir(parents=True, exist_ok=False)
    result = dict(status='complete_all_case_coordinate_archive_readback', terminal_state=state,
                  source_hashes=sources, checker_sha256=sha(__file__), entries=len(checked),
                  excluded_entries=excluded, archives=checked,
                  scope='All archive identities, compressed hashes, sizes and gzip CRCs verified; independent decompressed hashes recorded. Experimental methodology selection reconstructed from metadata. Atom-level residue correspondence and quality are not assessed.')
    (OUTPUT / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print('Complete:', len(checked), flush=True)


if __name__ == '__main__':
    main()
