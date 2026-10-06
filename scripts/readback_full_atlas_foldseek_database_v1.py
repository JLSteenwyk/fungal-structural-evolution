#!/usr/bin/env python3
"""Independently read the completed mixed-source Foldseek database."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np


def sha(path):
    """Hash arbitrarily large producer artifacts without whole-file allocation."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def parse_index(path):
    result = {}
    with Path(path).open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != 3:
                raise ValueError('Invalid native index row')
            key, offset, size = map(int, fields)
            if key in result or offset < 0 or size <= 0:
                raise ValueError('Repeated or invalid native index entry')
            result[key] = (offset, size)
    return result


def native_record(handle, location, text=False):
    offset, size = location
    handle.seek(offset)
    value = handle.read(size)
    if len(value) != size or not value.endswith(b'\0'):
        raise ValueError('Invalid native record boundary')
    payload = value[:-1]
    return payload.rstrip(b'\n') if text else payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer-receipt', required=True, type=Path)
    parser.add_argument('--producer-plan', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--receipt', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.receipt.exists():
        raise FileExistsError('Fresh readback paths required')
    plan = json.loads(args.producer_plan.read_text())
    producer = json.loads(args.producer_receipt.read_text())
    if producer['status'] != 'complete_full_prediction_atlas_foldseek_database_pending_independent_execution_and_database_readback':
        raise ValueError('Completed producer receipt required')
    if producer['plan_sha256'] != sha(args.producer_plan):
        raise ValueError('Producer plan digest differs')
    root = Path(plan['output']); prefix = root / 'structures'
    bindings = json.loads(Path(producer['source_and_output_bindings']).read_text())
    if sha(producer['source_and_output_bindings']) != producer['source_and_output_bindings_sha256']:
        raise ValueError('Changed producer binding map')
    for path, digest in bindings.items():
        if sha(path) != digest:
            raise ValueError('Changed bound producer input or output')
    expected = {}
    with Path(plan['models']).open() as handle:
        for line in handle:
            row = json.loads(line); model = row['model']; name = Path(model['path']).stem
            if name in expected or row['source'] not in ('AFDB', 'ESMFold'):
                raise ValueError('Repeated or unknown source model')
            expected[name] = (row['source'], model)
    lookup = {}
    with Path(str(prefix)+'.lookup').open() as handle:
        for line in handle:
            key, name, _ = line.rstrip('\n').split('\t')
            key = int(key)
            if key in lookup or name not in expected:
                raise ValueError('Foreign or repeated native alias')
            lookup[key] = name
    aa, ss, ca = (parse_index(str(prefix)+suffix) for suffix in ('.index', '_ss.index', '_ca.index'))
    if set(lookup) != set(aa) or set(lookup) != set(ss) or set(lookup) != set(ca) or len(lookup) != len(expected):
        raise ValueError('Native key grid differs')
    counts = {source: Counter(models=0, residues=0) for source in ('AFDB', 'ESMFold')}
    args.output.mkdir(parents=True)
    with prefix.open('rb') as amino, Path(str(prefix)+'_ss').open('rb') as states, Path(str(prefix)+'_ca').open('rb') as coords:
        for i, (key, name) in enumerate(lookup.items(), 1):
            source, model = expected[name]; sequence = native_record(amino, aa[key], text=True).upper()
            alphabet = native_record(states, ss[key], text=True).upper(); coordinate = native_record(coords, ca[key])
            length = model['length']
            if len(sequence) != length or hashlib.sha256(sequence).hexdigest() != model['sequence_sha256']:
                raise ValueError('Native amino-acid identity differs')
            if len(alphabet) != length or not set(alphabet) <= set(b'ACDEFGHIKLMNPQRSTVWYX'):
                raise ValueError('Native structural alphabet differs')
            if len(coordinate) != 12 * length or not np.isfinite(np.frombuffer(coordinate, dtype=np.float32)).all():
                raise ValueError('Native coordinate record differs')
            counts[source]['models'] += 1; counts[source]['residues'] += length
            if i % 100000 == 0:
                (args.output/'state.json').write_text(json.dumps({'stage':'independent_native_record_readback','models':i}, indent=2)+'\n')
    observed = {source: dict(value) for source, value in counts.items()}
    for source in observed:
        if observed[source]['models'] != plan['expected_models_by_source'][source] or observed[source]['residues'] != plan['expected_residues_by_source'][source]:
            raise ValueError('Source count differs')
    result = dict(status='passed_independent_full_prediction_atlas_foldseek_native_record_readback',
                  checked_utc=datetime.now(timezone.utc).isoformat(), counts_by_source=observed,
                  producer_receipt_sha256=sha(args.producer_receipt), producer_plan_sha256=sha(args.producer_plan),
                  independent_parser_sha256=sha(Path(__file__)), scientific_eligibility=False, gpu=False,
                  scope='Independent native index/record parser checks every full-atlas lookup identity, amino-acid hash, 3Di record length/alphabet, and float32 coordinate dimensions/finiteness. Original-CIF geometry and 3Di reconstruction remain producer-validated dependencies; no clustering, homology, confidence calibration or evolutionary inference.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    (args.output/'state.json').write_text(json.dumps({'stage':'complete','counts_by_source':observed}, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
