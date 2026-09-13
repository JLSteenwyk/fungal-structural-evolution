#!/usr/bin/env python3
"""Retrieve and validate the authors' published 3Di exchangeability models."""
import argparse
import hashlib
import json
import math
from pathlib import Path
from urllib.request import urlopen
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
MODELS = {'Q.3Di.AF': (311466, '8f897fffd14f247afbee44dad2fce2ce'),
          'Q.3Di.LLM': (311467, '3a6c99bdf5cc5fc72ba5eaedfc2dd2bd')}
API = 'https://edmond.mpg.de/api/datasets/:persistentId/?persistentId=doi:10.17617/3.1MJJBH'


def validate_model(raw):
    rows = [[float(x) for x in line.split()] for line in raw.decode('ascii').splitlines() if line.strip()]
    if [len(row) for row in rows] != list(range(1, 20)) + [20]:
        raise ValueError('Expected 190 lower-triangular exchangeabilities and 20 frequencies')
    if not all(math.isfinite(v) and v > 0 for row in rows for v in row):
        raise ValueError('Invalid exchangeability or frequency')
    # Frequencies were published rounded to six decimal places.
    if abs(sum(rows[-1]) - 1) > 20 * .5e-6 + 1e-12:
        raise ValueError('Frequencies do not sum to unity within publication rounding')
    return {'exchangeabilities': 190, 'frequencies': rows[-1], 'frequency_sum': sum(rows[-1])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable model directory')
    with urlopen(API, timeout=60) as response:
        deposit = response.read()
    version = json.loads(deposit)['data']['latestVersion']
    if (version['versionNumber'], version['versionMinorNumber']) != (3, 0):
        raise ValueError('Deposit version changed; review before repinning')
    entries = {f['dataFile']['id']: f['dataFile'] for f in version['files']}
    blobs, records = {}, {}
    for name, (file_id, expected_md5) in MODELS.items():
        entry = entries[file_id]
        if entry['filename'] != name or entry['checksum'] != {'type': 'MD5', 'value': expected_md5}:
            raise ValueError('Deposit identity/checksum differs')
        url = f'https://edmond.mpg.de/api/access/datafile/{file_id}'
        with urlopen(url, timeout=60) as response:
            raw = response.read()
        if len(raw) != entry['filesize'] or hashlib.md5(raw).hexdigest() != expected_md5:
            raise ValueError('Published file checksum or size mismatch')
        records[name] = validate_model(raw) | {'url': url, 'file_id': file_id,
            'publisher_md5': expected_md5, 'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw)}
        blobs[name] = raw
    args.output.mkdir(parents=True)
    blobs['deposit-v3.0.json'] = deposit
    for name, raw in blobs.items():
        (args.output / name).write_bytes(raw)
    receipt = {'status': 'complete_published_3di_model_validation',
        'doi': '10.17617/3.1MJJBH', 'paper_doi': '10.1093/molbev/msaf124',
        'dataset_version': '3.0', 'retrieved_utc': datetime.now(timezone.utc).isoformat(),
        'state_order': 'ARNDCQEGHILKMFPSTWYV',
        'format_reference': 'https://iqtree.github.io/doc/Substitution-Models#user-defined-empirical-protein-models',
        'interpretation': 'PAML lower-triangular symmetric exchangeabilities plus stationary frequencies; not an instantaneous rate matrix. Letters encode 3Di states, not amino acids. AF is the coordinate-derived primary model; LLM is a training-source sensitivity model.',
        'models': records, 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'artifacts': {name: hashlib.sha256(raw).hexdigest() for name, raw in blobs.items()}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
