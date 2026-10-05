#!/usr/bin/env python3
"""Independently decode every available original PAE matrix and reconstruct all summaries."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import multiprocessing
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def decode(entry):
    model, source = entry['model'], entry['source']
    if source == 'AFDB':
        cached = entry['cache']; receipt = cached['receipt']
        assert sha(cached['receipt_path']) == cached['receipt_sha256']
        assert json.loads(Path(cached['receipt_path']).read_text()) == receipt
        compressed = Path(receipt['path']).read_bytes()
        assert hashlib.sha256(compressed).hexdigest() == receipt['gzip_sha256']
        raw = gzip.decompress(compressed)
        assert len(raw) == receipt['json_bytes'] and hashlib.sha256(raw).hexdigest() == receipt['json_sha256']
        payload = json.loads(raw)
        assert isinstance(payload, list) and len(payload) == 1
        values = np.array(payload[0]['predicted_aligned_error'], dtype=np.float64)
        maximum = float(payload[0]['max_predicted_aligned_error'])
        assert math.isfinite(maximum) and maximum >= 0
        allowance = .51  # Original AFDB export-rounding contract, unchanged.
        size = len(compressed)
    else:
        assert source == 'ESMFold'
        receipt_path = Path(model['prediction_receipt_path']); matrix_path = Path(model['local_pae_npz_path'])
        assert sha(receipt_path) == model['prediction_receipt_sha256']
        assert sha(matrix_path) == model['local_pae_npz_sha256']
        receipt = json.loads(receipt_path.read_text())
        assert receipt['status'] == 'verified_prediction'
        assert receipt['sequence_sha256'] == model['sequence_sha256'] and receipt['length'] == model['length']
        assert receipt['config_sha256'] == model['prediction_config_sha256']
        assert receipt['artifacts'][matrix_path.name] == model['local_pae_npz_sha256']
        with np.load(matrix_path, allow_pickle=False) as payload:
            sequence = str(payload['sequence']); values = payload['pae'].copy()
            maximum = float(payload['max_pae'])
        assert len(sequence) == model['length'] and hashlib.sha256(sequence.encode()).hexdigest() == model['sequence_sha256']
        assert math.isfinite(maximum) and maximum > 0 and maximum == receipt['max_predicted_aligned_error']
        allowance = 1e-4  # Original ESMFold exporter contract, unchanged.
        size = matrix_path.stat().st_size
    assert values.shape == (model['length'], model['length'])
    assert np.isfinite(values).all() and np.min(values) >= 0 and np.max(values) <= maximum+allowance
    return values, maximum, size


def reconstruct(entry, observed, values, declared, size):
    model, source = entry['model'], entry['source']
    # Sort-based order statistics and extended precision reduction are distinct
    # from the producer's float64 sum and boolean-array threshold reductions.
    ordered = np.sort(values.reshape(-1))
    mean = float(np.sum(values, dtype=np.longdouble)/np.longdouble(values.size))
    expected = dict(source=source, model_id=model['model_id'], version=model['version'],
                    sequence_sha256=model['sequence_sha256'], length=model['length'],
                    directional_entries=int(values.size), source_bytes=size, original_dtype=values.dtype.str,
                    numeric_sha256=hashlib.sha256(values.dtype.str.encode()+str(values.shape).encode()+values.tobytes()).hexdigest(),
                    minimum=float(ordered[0]), maximum=float(ordered[-1]), declared_maximum=declared,
                    mean=mean, pae_le5=int(np.searchsorted(ordered, 5, side='right')),
                    pae_le10=int(np.searchsorted(ordered, 10, side='right')),
                    pae_le20=int(np.searchsorted(ordered, 20, side='right')),
                    status='validated_original_directional_matrix', symmetrized=False)
    assert set(observed) == set(expected)
    for key, value in expected.items():
        if key != 'mean': assert observed[key] == value, key
    assert math.isfinite(observed['mean'])
    # A declared roundoff bound for nonnegative pairwise float64 summation,
    # separate from unchanged PAE bounds and all evolutionary fit tolerances.
    bound = np.finfo(np.float64).eps*(8*math.ceil(math.log2(max(2, values.size)))+16)*max(1.0, abs(mean))
    error = abs(observed['mean']-mean)
    assert error <= bound, (error, bound)
    return error, bound


def read_job(proof):
    assert sha(proof['job']) == proof['job_sha256'] and sha(proof['output']) == proof['output_sha256']
    counts, cells, sizes = Counter(), Counter(), Counter()
    max_error = 0.0; max_ratio = 0.0
    with Path(proof['job']).open() as inputs, Path(proof['output']).open() as reports:
        rows = (json.loads(line) for line in reports)
        for line in inputs:
            entry = json.loads(line); observed = next(rows)
            matrix, maximum, size = decode(entry)
            error, bound = reconstruct(entry, observed, matrix, maximum, size)
            source = entry['source']; counts[source] += 1; cells[source] += matrix.size; sizes[source] += size
            max_error = max(max_error, error); max_ratio = max(max_ratio, error/bound)
        assert next(rows, None) is None
    assert dict(counts) == proof['models'] and dict(cells) == proof['directional_entries']
    assert dict(sizes) == proof['source_bytes']
    return dict(job=proof['job'], matrices=dict(counts), directional_entries=dict(cells), source_bytes=dict(sizes),
                maximum_absolute_mean_error=max_error, maximum_mean_roundoff_bound_fraction=max_ratio)


def self_test():
    values = np.array([[0.25, 4], [12.5, 0.25]], dtype=np.float64)
    entry = dict(source='AFDB', model=dict(model_id='OFFLINE', version=6, sequence_sha256='literal', length=2))
    observed = dict(source='AFDB', model_id='OFFLINE', version=6, sequence_sha256='literal', length=2,
                    directional_entries=4, source_bytes=0, original_dtype='<f8',
                    numeric_sha256=hashlib.sha256(b'<f8'+b'(2, 2)'+values.tobytes()).hexdigest(),
                    minimum=.25, maximum=12.5, declared_maximum=12.5, mean=4.25,
                    pae_le5=3, pae_le10=3, pae_le20=4, status='validated_original_directional_matrix', symmetrized=False)
    assert reconstruct(entry, observed, values, 12.5, 0)[0] == 0
    checks = ['literal_directional_values_and_scalar_statistics']
    for key, value in [('pae_le5', 4), ('mean', 4.251), ('numeric_sha256', 'changed'),
                       ('symmetrized', True), ('original_dtype', '<f4')]:
        try: reconstruct(entry, dict(observed, **{key:value}), values, 12.5, 0)
        except AssertionError: checks.append('reject_'+key)
        else: raise AssertionError('Corruption accepted: '+key)
    assert np.finfo(np.longdouble).eps < np.finfo(np.float64).eps
    checks.append('extended_precision_available_on_runtime')
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args(); assert not args.receipt.exists()
    checks = self_test()
    if args.self_test:
        result = dict(status='passed_literal_independent_pae_statistic_controls', checks=checks,
                      source_hashes={str(Path(__file__)):sha(Path(__file__))}, scientific_eligibility=False)
    else:
        assert args.plan is not None
        plan = json.loads(args.plan.read_text()); pins = dict(plan['pins']); verify(pins)
        pp, ep = Path(plan['producer_receipt']), Path(plan['recovered_native_evidence'])
        producer, evidence = [json.loads(p.read_text()) for p in (pp, ep)]
        assert producer['status'] == 'passed_full_atlas_pae_availability_replay_and_all_available_matrix_validation'
        assert evidence['status'] == 'verified_original_native_receipt_and_whole_journal_without_api_wait'
        assert evidence['native_exit_code'] == 0 and evidence['validation_sha256'] == sha(pp)
        assert evidence['wrapper_initial_and_terminal_matched']
        assert evidence['manager_start_records'] == evidence['manager_completion_records'] == 1
        assert evidence['original_tool_terminal_available'] is False and evidence['original_tool_terminal_exit_code'] is None
        root = Path(plan['output']); assert not root.exists(); root.mkdir(parents=True)
        assert len(producer['proofs']) == producer['jobs'] and len({p['job'] for p in producer['proofs']}) == producer['jobs']
        counts, cells, sizes = Counter(), Counter(), Counter(); proofs = []
        with ProcessPoolExecutor(max_workers=plan['cpu'], mp_context=multiprocessing.get_context('spawn')) as pool:
            futures = [pool.submit(read_job, proof) for proof in producer['proofs']]
            for future in as_completed(futures):
                proof = future.result(); proofs.append(proof)
                counts.update(proof['matrices']); cells.update(proof['directional_entries']); sizes.update(proof['source_bytes'])
                tmp = root/'state.tmp'; tmp.write_text(json.dumps(dict(completed_jobs=len(proofs), total_jobs=producer['jobs'], matrices=dict(counts)), indent=2)+'\n'); tmp.replace(root/'state.json')
        assert dict(counts) == producer['matrices_validated'] and dict(cells) == producer['directional_entries']
        assert dict(sizes) == producer['source_matrix_bytes']
        for mapping in (producer['source_hashes'], evidence['source_hashes']):
            for path, digest in mapping.items(): bind(pins, path, digest)
        for path in (args.plan, pp, ep, Path(__file__)): bind(pins, path)
        verify(pins)
        result = dict(status='passed_full_independent_available_pae_decode_and_statistic_reconstruction',
                      checked_utc=datetime.now(timezone.utc).isoformat(), checks=checks, matrices=dict(counts),
                      directional_entries=dict(cells), source_matrix_bytes=dict(sizes),
                      jobs=len(proofs), proofs=sorted(proofs, key=lambda p:p['job']),
                      maximum_absolute_mean_error=max(p['maximum_absolute_mean_error'] for p in proofs),
                      maximum_mean_roundoff_bound_fraction=max(p['maximum_mean_roundoff_bound_fraction'] for p in proofs),
                      original_producer_api_terminal_available=False, source_hashes=pins,
                      scientific_eligibility=False, gpu=False, new_predictions=0,
                      scope='Every available original AFDB/ESMFold matrix decoded without producer decoder imports; '
                            'all directional numeric bytes, identities, extrema and threshold counts reconstructed. '
                            'Extended precision means use a declared positive-sum roundoff bound. No symmetrization, '
                            'new PAE value bounds, prediction, calibrated accuracy, context boundaries or biological '
                            'effects. All missing matrices remain explicit. Original producer native completion is '
                            'journal-proven; lost original API terminal remains unavailable.')
    with args.receipt.open('x') as handle: json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes', 'proofs')}, indent=2))


if __name__ == '__main__': main()
