#!/usr/bin/env python3
"""Freeze the entire PAE cache and enumerate exact-version availability for every atlas model."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import fcntl
import gzip
import io
import json
import os
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def classify(raw, cache):
    source, model = raw['source'], raw['model']
    if source == 'ESMFold':
        return 'local_prediction_candidate_pending_matrix_validation', None
    assert source == 'AFDB'
    key = (model['model_id'], model['version'])
    cached = cache.get(key)
    expected_url = f"https://alphafold.ebi.ac.uk/files/{model['model_id']}-predicted_aligned_error_v{model['version']}.json"
    if cached is not None:
        r = cached['receipt']
        matches = all(r[k] == model[k] for k in ('model_id', 'version', 'length', 'sequence_sha256'))
        matches = matches and r['url'] == expected_url and model.get('pae_url') == expected_url
        return ('exact_cached_candidate_pending_matrix_validation' if matches else 'cache_provenance_mismatch'), cached
    return ('missing_cached_matrix_retrievable' if model.get('pae_url') == expected_url
            else 'missing_cached_matrix_without_matching_advertised_url'), None


def compressed_text(path):
    raw = path.open('xb')
    gz = gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0)
    # Closing the returned text stream closes gzip; the original descriptor is
    # explicitly owned by the caller and closed separately.
    return raw, io.TextIOWrapper(gz, encoding='utf-8', newline='')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args(); assert not a.receipt.exists()
    plan = json.loads(a.plan.read_text()); pins = dict(plan['pins']); verify(pins)
    atlas = json.loads(Path(plan['atlas_closure']).read_text())
    assert atlas['status'] == 'complete_verified_full_prediction_atlas_availability_union'
    assert atlas['model_counts'] == plan['expected_models_by_source']
    cache_root = Path(plan['cache']); output = Path(plan['output']); assert not output.exists()
    output.mkdir(parents=True)
    cache, matched, counts, sources = {}, set(), Counter(), Counter()
    matrix_cells = Counter(); missing_cells = 0; cache_bytes = 0
    lock = (cache_root/'.retrieval.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    try:
        names = {e.name for e in os.scandir(cache_root) if e.is_file()}
        matrices = {name for name in names if name.endswith('.json.gz')}
        receipts = {name for name in names if name.endswith('.receipt.json')}
        assert matrices == {name[:-len('.receipt.json')]+'.json.gz' for name in receipts}
        assert names-matrices-receipts == {'.retrieval.lock'}
        cache_path = output/'cache_metadata.jsonl'
        with cache_path.open('x') as emitted:
            for name in sorted(receipts):
                path = cache_root/name; r = json.loads(path.read_text())
                key = (r['model_id'], r['version']); assert key not in cache
                assert r['status'] == 'verified'
                label = f"{r['model_id']}-v{r['version']}"
                assert name == label+'.receipt.json'
                matrix = cache_root/(label+'.json.gz')
                assert Path(r['path']).resolve() == matrix.resolve()
                assert matrix.stat().st_size == r['compressed_bytes']
                bind(pins, matrix, r['gzip_sha256']); bind(pins, path)
                record = dict(receipt_path=str(path), receipt_sha256=sha(path), receipt=r)
                cache[key] = record; cache_bytes += r['compressed_bytes']
                emitted.write(json.dumps(record, separators=(',', ':'), allow_nan=False)+'\n')
        # Hash all cached files under the original retrieval mutex, including
        # entries that do not match any currently selected source model.
        verify(pins)
        full_path, queue_path = output/'model_pae_dispositions.jsonl.gz', output/'missing_afdb_models.jsonl'
        file_handle, emitted = compressed_text(full_path)
        try:
            with emitted, queue_path.open('x') as queue, Path(plan['models']).open() as models:
                for line in models:
                    raw = json.loads(line); source, m = raw['source'], raw['model']
                    status, record = classify(raw, cache)
                    sources[source] += 1; counts[source+':'+status] += 1
                    entry = dict(source=source, model=m, status=status)
                    if record is not None:
                        entry['cache'] = record; matched.add((m['model_id'], m['version']))
                    if status == 'missing_cached_matrix_retrievable':
                        queue.write(json.dumps(m, separators=(',', ':'), allow_nan=False)+'\n')
                        missing_cells += m['length']**2
                    if status in ('exact_cached_candidate_pending_matrix_validation', 'local_prediction_candidate_pending_matrix_validation'):
                        matrix_cells[source] += m['length']**2
                    emitted.write(json.dumps(entry, separators=(',', ':'), allow_nan=False)+'\n')
        finally:
            file_handle.close()
        assert dict(sources) == plan['expected_models_by_source']
        final_names = {e.name for e in os.scandir(cache_root) if e.is_file()}
        assert final_names == names
        verify(pins)
    finally:
        fcntl.flock(lock, fcntl.LOCK_UN); lock.close()
    raw_result = dict(status='complete_full_atlas_pae_cache_and_exact_model_availability_snapshot',
                      sources=dict(sources), dispositions=dict(counts), cache_entries=len(cache),
                      cache_matrix_bytes=cache_bytes, cache_entries_matching_any_selected_model=len(matched),
                      cache_entries_outside_current_model_inventory=len(cache)-len(matched),
                      candidate_directional_matrix_cells=dict(matrix_cells), missing_retrievable_matrix_cells=missing_cells,
                      plan_sha256=sha(a.plan), artifacts={p.name: sha(p) for p in (cache_path, full_path, queue_path)},
                      scope='Entire original AFDB cache receipt/matrix census and byte bindings under original '
                            'exclusive retrieval mutex; every atlas source model classified by exact model/version/sequence/length/URL. '
                            'All local ESMFold candidates and missing/unmatched/conflicting records retained. '
                            'Matrix content, full independent availability readback, confidence calibration '
                            'and evolutionary inference are not qualified by this snapshot.')
    raw_path = output/'receipt.json'
    with raw_path.open('x') as handle:
        json.dump(raw_result, handle, indent=2, allow_nan=False); handle.write('\n')
    for path in (a.plan, raw_path, cache_path, full_path, queue_path, Path(__file__)):
        bind(pins, path)
    verify(pins)
    result = dict(raw_result, checked_utc=datetime.now(timezone.utc).isoformat(), source_hashes=pins,
                  raw_receipt_sha256=sha(raw_path), scientific_eligibility=False, gpu=False, new_predictions=0)
    with a.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_hashes', 'artifacts')}, indent=2))


if __name__ == '__main__':
    main()
