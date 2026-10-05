#!/usr/bin/env python3
"""Replay every full-atlas PAE availability row and validate all available original matrices."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import json
import multiprocessing
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify
from retrieve_marker_pae import validate_pae


def load_matrix(entry):
    source, model = entry['source'], entry['model']
    if source == 'AFDB':
        cached = entry['cache']; r = cached['receipt']
        rp = Path(cached['receipt_path']); assert sha(rp) == cached['receipt_sha256']
        assert json.loads(rp.read_text()) == r
        compressed = Path(r['path']).read_bytes()
        assert hashlib.sha256(compressed).hexdigest() == r['gzip_sha256']
        raw = gzip.decompress(compressed)
        assert hashlib.sha256(raw).hexdigest() == r['json_sha256'] and len(raw) == r['json_bytes']
        # Qualified monomer parser retains the original +0.51 export-rounding bound.
        matrix = validate_pae(raw, model['length'])
        maximum = float(json.loads(raw)[0]['max_predicted_aligned_error'])
        source_bytes = len(compressed)
    else:
        assert source == 'ESMFold'
        pr = Path(model['prediction_receipt_path']); npz = Path(model['local_pae_npz_path'])
        assert sha(pr) == model['prediction_receipt_sha256'] and sha(npz) == model['local_pae_npz_sha256']
        predicted = json.loads(pr.read_text())
        assert predicted['status'] == 'verified_prediction'
        assert predicted['sequence_sha256'] == model['sequence_sha256'] and predicted['length'] == model['length']
        assert predicted['config_sha256'] == model['prediction_config_sha256']
        assert predicted['artifacts'][npz.name] == model['local_pae_npz_sha256']
        with np.load(npz, allow_pickle=False) as values:
            sequence = str(values['sequence']); matrix = values['pae'].copy(); maximum = float(values['max_pae'])
        assert len(sequence) == model['length'] and hashlib.sha256(sequence.encode()).hexdigest() == model['sequence_sha256']
        assert matrix.shape == (model['length'], model['length'])
        assert np.isfinite(matrix).all() and not (matrix < 0).any()
        assert np.isfinite(maximum) and maximum > 0 and not (matrix > maximum+1e-4).any()
        assert maximum == predicted['max_predicted_aligned_error']
        source_bytes = npz.stat().st_size
    return matrix, maximum, source_bytes


def audit_job(job, output):
    job = Path(job); target = Path(output)/(job.stem+'.jsonl')
    assert not target.exists()
    counts = Counter(); cells = Counter(); source_bytes = Counter()
    with job.open() as inp, target.open('x') as out:
        for line in inp:
            entry = json.loads(line); source, model = entry['source'], entry['model']
            matrix, maximum, size = load_matrix(entry)
            row = dict(source=source, model_id=model['model_id'], version=model['version'],
                       sequence_sha256=model['sequence_sha256'], length=model['length'],
                       directional_entries=int(matrix.size), source_bytes=size, original_dtype=matrix.dtype.str,
                       numeric_sha256=hashlib.sha256(matrix.dtype.str.encode()+str(matrix.shape).encode()+matrix.tobytes()).hexdigest(),
                       minimum=float(matrix.min()), maximum=float(matrix.max()), declared_maximum=maximum,
                       mean=float(np.sum(matrix, dtype=np.float64)/matrix.size),
                       pae_le5=int((matrix <= 5).sum()), pae_le10=int((matrix <= 10).sum()),
                       pae_le20=int((matrix <= 20).sum()),
                       status='validated_original_directional_matrix', symmetrized=False)
            counts[source] += 1; cells[source] += matrix.size; source_bytes[source] += size
            out.write(json.dumps(row, sort_keys=True, allow_nan=False)+'\n')
    return dict(job=str(job), job_sha256=sha(job), output=str(target), output_sha256=sha(target),
                models=dict(counts), directional_entries=dict(cells), source_bytes=dict(source_bytes))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text()); pins = dict(plan['pins']); verify(pins)
    producer_path, tp = Path(plan['producer_receipt']), Path(plan['producer_transport'])
    producer, transport = [json.loads(p.read_text()) for p in (producer_path, tp)]
    assert producer['status'] == 'complete_full_atlas_pae_cache_and_exact_model_availability_snapshot'
    assert transport['status'] == 'verified_original_software_wait_and_whole_wrapper_payloads'
    assert transport['validation_sha256'] == sha(producer_path) and transport['original_tool_terminal_exit_code'] == 0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    assert transport['manager_start_records'] == transport['manager_completion_records'] == 1
    for mapping in (producer['source_hashes'], transport['source_hashes']):
        for path, digest in mapping.items():bind(pins, path, digest)
    verify(pins)
    root = Path(plan['inventory']); output = Path(plan['output']); assert not output.exists()
    output.mkdir(parents=True); jobs=output/'jobs'; reports=output/'reports'; jobs.mkdir(); reports.mkdir()
    cache = {}
    with (root/'cache_metadata.jsonl').open() as inp:
        for line in inp:
            entry = json.loads(line); r = entry['receipt']; key=(r['model_id'], r['version'])
            assert key not in cache and sha(entry['receipt_path']) == entry['receipt_sha256']
            assert json.loads(Path(entry['receipt_path']).read_text()) == r
            assert sha(r['path']) == r['gzip_sha256']; cache[key]=entry
    assert len(cache) == producer['cache_entries']
    missing, counts, sources = 0, Counter(), Counter(); buffer=[]; job_paths=[]
    def flush():
        path=jobs/f'shard_{len(job_paths):05d}.jsonl'
        with path.open('x') as handle:
            for row in buffer:handle.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
        job_paths.append(path);buffer.clear()
    with gzip.open(root/'model_pae_dispositions.jsonl.gz', 'rt') as emitted, Path(plan['models']).open() as original, (root/'missing_afdb_models.jsonl').open() as queue:
        rows=(json.loads(line) for line in emitted); requested=(json.loads(line) for line in queue)
        for line in original:
            source_row=json.loads(line);source,m=source_row['source'],source_row['model'];row=next(rows)
            assert row['source']==source and row['model']==m
            available=False
            if source=='ESMFold':
                expected='local_prediction_candidate_pending_matrix_validation';available=True
                assert 'cache' not in row
            else:
                assert source=='AFDB';record=cache.get((m['model_id'],m['version']))
                url=f"https://alphafold.ebi.ac.uk/files/{m['model_id']}-predicted_aligned_error_v{m['version']}.json"
                if record is None:
                    assert 'cache' not in row
                    expected=('missing_cached_matrix_retrievable' if m.get('pae_url')==url
                              else 'missing_cached_matrix_without_matching_advertised_url')
                    if expected=='missing_cached_matrix_retrievable':assert next(requested)==m;missing+=1
                else:
                    assert row['cache']==record
                    identity=(record['receipt']['sequence_sha256'],record['receipt']['length'],record['receipt']['url'])
                    available=identity==(m['sequence_sha256'],m['length'],url) and m.get('pae_url')==url
                    expected=('exact_cached_candidate_pending_matrix_validation' if available else 'cache_provenance_mismatch')
            assert row['status']==expected
            sources[source]+=1;counts[source+':'+expected]+=1
            if available:
                buffer.append(row)
                if len(buffer)==plan['models_per_job']:flush()
        assert next(rows,None) is None and next(requested,None) is None
    if buffer:flush()
    assert dict(sources)==producer['sources'] and dict(counts)==producer['dispositions']
    totals, cells, sizes, proofs = Counter(), Counter(), Counter(), []
    with ProcessPoolExecutor(max_workers=plan['resources']['cpu'], mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[pool.submit(audit_job,str(job),str(reports)) for job in job_paths]
        for future in as_completed(futures):
            try:r=future.result()
            except BaseException:
                for pending in futures:pending.cancel()
                raise
            proofs.append(r);totals.update(r['models']);cells.update(r['directional_entries']);sizes.update(r['source_bytes'])
            temporary=output/'state.tmp'
            temporary.write_text(json.dumps(dict(stage='validating_available_matrices',completed_jobs=len(proofs),total_jobs=len(job_paths),models=dict(totals),directional_entries=dict(cells)),indent=2)+'\n')
            temporary.replace(output/'state.json')
    expected_models=dict(AFDB=counts['AFDB:exact_cached_candidate_pending_matrix_validation'],ESMFold=counts['ESMFold:local_prediction_candidate_pending_matrix_validation'])
    assert dict(totals)==expected_models and dict(cells)==producer['candidate_directional_matrix_cells']
    for proof in proofs:
        bind(pins,proof['job'],proof['job_sha256']);bind(pins,proof['output'],proof['output_sha256'])
    for path in (args.plan,producer_path,tp,Path(__file__)):bind(pins,path)
    verify(pins)
    result=dict(status='passed_full_atlas_pae_availability_replay_and_all_available_matrix_validation',
                checked_utc=datetime.now(timezone.utc).isoformat(),source_models=dict(sources),dispositions=dict(counts),
                matrices_validated=dict(totals),directional_entries=dict(cells),source_matrix_bytes=dict(sizes),
                missing_retrievable_afdb_models=missing,jobs=len(proofs),proofs=sorted(proofs,key=lambda p:p['job']),
                source_hashes=pins,scientific_eligibility=False,gpu=False,new_predictions=0,
                scope='Every full source-model disposition and missing queue row independently replayed; '
                      'all available exact-version AFDB and original ESMFold matrices byte/provenance/shape/finite/max-bound '
                      'validated with all directional values retained and hashed, no symmetrization. '
                      'AFDB monomer validator shared with original retrieval. Missing matrices, independent '
                      'matrix-statistic reconstruction, native context confidence, accuracy calibration and '
                      'biological inference remain separate requirements.')
    with args.receipt.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes','proofs')},indent=2))


if __name__=='__main__':main()
