#!/usr/bin/env python3
"""Audit every source model and preserve original C-alpha coordinates/confidence in tar shards."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import io
import json
import multiprocessing
from pathlib import Path
import shutil
import sys
import tarfile
import time

import Bio
import numpy as np

from ancestral_chain_attempt import sha
from extract_domain_coordinates import load_atoms
from reference_measurement_union_sources import bind, verify


COUNT_KEYS = ('models', 'valid_models', 'rejected_models', 'source_coordinate_bytes',
              'valid_residues', 'missing_backbone_residues', 'plddt_ge50', 'plddt_ge70', 'plddt_ge90')


def array_digest(value):
    return hashlib.sha256(value.dtype.str.encode()+str(value.shape).encode()+value.tobytes()).hexdigest()


def profile(model):
    sequence, residues, confidence = load_atoms(model)
    n = len(sequence)
    arrays = dict(sequence_ascii=np.frombuffer(sequence.encode('ascii'), dtype=np.uint8).copy(),
                  ca_xyz=np.asarray([next(a[2:5] for a in residues[p] if a[0] == 'CA')
                                     for p in range(1, n+1)], dtype=np.float64),
                  ca_plddt=np.asarray([confidence[p] for p in range(1, n+1)], dtype=np.float64),
                  missing_backbone=np.asarray([not {'N', 'CA', 'C', 'O'} <= {a[0] for a in residues[p]}
                                                for p in range(1, n+1)], dtype=np.bool_))
    assert arrays['ca_xyz'].shape == (n, 3) and arrays['ca_plddt'].shape == (n,)
    assert np.isfinite(arrays['ca_xyz']).all() and np.isfinite(arrays['ca_plddt']).all()
    summary = dict(residues=n, missing_backbone_residues=int(arrays['missing_backbone'].sum()),
                   mean_coordinate_ca_plddt=sum(confidence.values())/n,
                   minimum_ca_plddt=float(arrays['ca_plddt'].min()),
                   maximum_ca_plddt=float(arrays['ca_plddt'].max()))
    for cutoff in (50, 70, 90):
        summary['plddt_ge'+str(cutoff)] = int((arrays['ca_plddt'] >= cutoff).sum())
    summary['catalog_mean_ca_plddt'] = model['mean_ca_plddt']
    summary['coordinate_minus_catalog_mean_ca_plddt'] = summary['mean_coordinate_ca_plddt']-model['mean_ca_plddt']
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **arrays)
    blob = buffer.getvalue()
    with np.load(io.BytesIO(blob), allow_pickle=False) as readback:
        assert set(readback.files) == set(arrays)
        assert all(np.array_equal(readback[key], value) for key, value in arrays.items())
    summary['array_sha256'] = {key: array_digest(value) for key, value in arrays.items()}
    return blob, summary


def shard(job_path, out_dir, reserve_gib):
    job_path, out = Path(job_path), Path(out_dir)
    label = job_path.stem
    tar_path, rows_path, receipt_path = [out/(label+s) for s in ('.tar', '.jsonl', '.receipt.json')]
    assert not any(p.exists() for p in (tar_path, rows_path, receipt_path))
    job_hash = sha(job_path)
    counts = {source: Counter(dict.fromkeys(COUNT_KEYS, 0)) for source in ('AFDB', 'ESMFold')}
    with tarfile.open(tar_path, 'w', format=tarfile.USTAR_FORMAT) as archive, rows_path.open('x') as rows:
        with job_path.open() as handle:
            for line in handle:
                if shutil.disk_usage(out).free < reserve_gib*2**30:
                    raise RuntimeError('Emergency disk reserve reached')
                raw = json.loads(line)
                source, model = raw['source'], raw['model']
                c = counts[source]
                c['models'] += 1
                c['source_coordinate_bytes'] += Path(model['path']).stat().st_size
                row = dict(source=source, sequence_sha256=model['sequence_sha256'], model_id=model['model_id'],
                           version=model['version'], source_path=model['path'], source_sha256=model['sha256'],
                           length=model['length'], prediction_config_sha256=model.get('prediction_config_sha256', ''),
                           confidence_representation='Original mmCIF C-alpha B_iso_or_equiv; score, not calibrated accuracy',
                           pae_evaluated=False)
                try:
                    blob, summary = profile(model)
                except (ValueError, KeyError) as error:
                    # Immutable byte corruption is a stage failure, never an eligibility exclusion.
                    if str(error) == 'Source coordinate checksum mismatch':
                        raise
                    row.update(status='source_validation_rejected', reason=type(error).__name__+': '+str(error))
                    c['rejected_models'] += 1
                else:
                    member = source+'-'+model['sequence_sha256']+'.npz'
                    info = tarfile.TarInfo(member)
                    info.size, info.mtime, info.mode = len(blob), 0, 0o644
                    archive.addfile(info, io.BytesIO(blob))
                    row.update(status='coordinate_and_residue_profile_exported', member=member,
                               npz_sha256=hashlib.sha256(blob).hexdigest(), **summary)
                    c['valid_models'] += 1
                    c['valid_residues'] += summary['residues']
                    for key in ('missing_backbone_residues', 'plddt_ge50', 'plddt_ge70', 'plddt_ge90'):
                        c[key] += summary[key]
                rows.write(json.dumps(row, sort_keys=True, allow_nan=False)+'\n')
    result = dict(job_sha256=job_hash, counts={s: dict(c) for s, c in counts.items()},
                  artifacts={p.name: sha(p) for p in (tar_path, rows_path)})
    with receipt_path.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text())
    assert '.'.join(map(str, sys.version_info[:3])) == plan['environment']['python']
    assert Bio.__version__ == plan['environment']['biopython'] and np.__version__ == plan['environment']['numpy']
    pins = dict(plan['pins'])
    verify(pins)
    closure = json.loads(Path(plan['atlas_closure']).read_text())
    assert closure['status'] == 'complete_verified_full_prediction_atlas_availability_union'
    assert closure['model_counts'] == plan['expected_models_by_source']
    assert closure['representative_proteins'] == 5815847 and closure['taxa'] == 526
    assert all(t['actual_terminal_exit_code'] == 0 and t['whole_original_payloads_verified']
               for t in closure['original_transports'])
    root = Path(plan['output'])
    assert not root.exists()
    assert shutil.disk_usage(root.parent).free >= plan['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(parents=True)
    jobs, shards = root/'jobs', root/'shards'
    jobs.mkdir(); shards.mkdir()
    job_paths, buffer, sources = [], [], Counter()
    def write_job():
        target = jobs/f'shard_{len(job_paths):05d}.jsonl'
        with target.open('x') as handle:
            for row in buffer:
                handle.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
        job_paths.append(target)
    with Path(plan['models']).open() as handle:
        for line in handle:
            row = json.loads(line)
            assert row['source'] in plan['expected_models_by_source']
            sources[row['source']] += 1
            buffer.append(row)
            if len(buffer) == plan['models_per_shard']:
                write_job(); buffer.clear()
    if buffer:
        write_job(); buffer.clear()
    assert dict(sources) == plan['expected_models_by_source']
    counts = {source: Counter(dict.fromkeys(COUNT_KEYS, 0)) for source in sources}
    proofs, started = [], time.monotonic()
    def state(status):
        value = dict(stage=status, completed_shards=len(proofs), total_shards=len(job_paths),
                     counts={s: dict(c) for s, c in counts.items()}, elapsed_seconds=time.monotonic()-started)
        temporary = root/'state.tmp'
        temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
        temporary.replace(root/'state.json')
    state('audit_running')
    with ProcessPoolExecutor(max_workers=plan['resources']['cpu'], mp_context=multiprocessing.get_context('spawn')) as pool:
        futures = {pool.submit(shard, str(job), str(shards), plan['resources']['emergency_free_disk_gib']): job for job in job_paths}
        for future in as_completed(futures):
            try:
                result = future.result()
            except BaseException:
                for pending in futures:
                    pending.cancel()
                raise
            job = futures[future]
            proofs.append(dict(job=str(job), receipt=result))
            for source, local in result['counts'].items():
                counts[source].update(local)
            state('audit_running')
    for source, c in counts.items():
        assert c['models'] == sources[source] and c['valid_models']+c['rejected_models'] == sources[source]
        assert 0 <= c['plddt_ge90'] <= c['plddt_ge70'] <= c['plddt_ge50'] <= c['valid_residues']
    verify(pins)
    raw = dict(status='complete_full_atlas_coordinate_and_residue_dispositions_pending_independent_readback',
               counts={s: dict(c) for s, c in counts.items()}, shards=len(proofs),
               plan_sha256=sha(args.plan), proofs=sorted(proofs, key=lambda x: x['job']),
               scope='Every source model attempted with unchanged qualified atom checks. Original float64 '
                     'C-alpha coordinates/confidence and exact sequence exported without rounding; missing '
                     'backbone and rejection dispositions retained. Inclusive score50/70/90 counts are '
                     'descriptive sensitivity masks, not acceptance/calibration. PAE/context features, '
                     'independent archive/atom/rejection readback and biological inference remain required.')
    raw_path = root/'receipt.json'
    with raw_path.open('x') as handle:
        json.dump(raw, handle, indent=2, allow_nan=False); handle.write('\n')
    state('native_audit_complete_output_binding_running')
    for proof in raw['proofs']:
        bind(pins, proof['job'], proof['receipt']['job_sha256'])
        for name, digest in proof['receipt']['artifacts'].items():
            bind(pins, shards/name, digest)
        bind(pins, shards/(Path(proof['job']).stem+'.receipt.json'))
    for path in (args.plan, raw_path, Path(__file__)):
        bind(pins, path)
    verify(pins)
    result = dict(status='completed_full_atlas_coordinate_and_residue_audit_pending_independent_readback',
                  checked_utc=datetime.now(timezone.utc).isoformat(), counts=raw['counts'], shards=raw['shards'],
                  raw_receipt_sha256=sha(raw_path), source_hashes=pins,
                  scientific_eligibility=False, gpu=False, new_predictions=0, scope=raw['scope'])
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    state('producer_complete_pending_original_wait_and_independent_readback')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
