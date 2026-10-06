#!/usr/bin/env python3
"""Encode every current source model and compare every native C-alpha to its CIF."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import multiprocessing
import os
from pathlib import Path
import re
import subprocess
import time

import numpy as np
import psutil

from ancestral_chain_attempt import sha
from build_whole_proteome_foldseek_database import index, readback
from readback_full_atlas_coordinate_profiles_v1 import original_profile
from reference_measurement_union_sources import bind, verify


def model_record(raw):
    if set(raw) != {'source', 'model'} or raw['source'] not in ('AFDB', 'ESMFold'):
        raise ValueError('Unknown source or model-record schema')
    model = raw['model']
    if type(model['length']) is not int or model['length'] <= 0:
        raise ValueError('Invalid model length')
    for key in ('sha256', 'sequence_sha256'):
        if not isinstance(model[key], str) or not re.fullmatch('[0-9a-f]{64}', model[key]):
            raise ValueError('Invalid immutable model digest')
    path = Path(model['path'])
    if path.suffix != '.cif' or '\n' in str(path) or '\t' in str(path):
        raise ValueError('Invalid CIF path-list input')
    return dict(source=raw['source'], model={key: model[key] for key in
        ('path', 'sha256', 'sequence_sha256', 'length', 'model_id', 'version')})


def unique_alias(row, aliases):
    name = Path(row['model']['path']).stem
    if not name or name in aliases:
        raise ValueError('Ambiguous source-model filename')
    aliases.add(name)
    return name


def compare_coordinates(raw, model):
    expected = original_profile(model)
    n = model['length']
    if len(raw) != 12 * n + 1 or raw[-1:] != b'\0':
        raise ValueError('Invalid native coordinate record')
    actual = np.frombuffer(raw[:-1], dtype=np.float32).reshape(3, n).T
    source = expected['ca_xyz'].astype(np.float32)
    if not np.isfinite(actual).all() or not np.array_equal(actual, source):
        raise ValueError('Native C-alpha coordinates differ from original CIF')
    return n, int(expected['missing_backbone'].sum())


def geometry_shard(job_path, job_hash, coordinate_path, output):
    job_path, output = Path(job_path), Path(output)
    if sha(job_path) != job_hash:
        raise ValueError('Geometry job changed')
    counts = {source: dict(models=0, residues=0, source_cif_bytes=0,
                          missing_backbone_residues=0) for source in ('AFDB', 'ESMFold')}
    with Path(coordinate_path).open('rb') as coordinates, job_path.open() as jobs:
        for line in jobs:
            row = json.loads(line)
            coordinates.seek(row['offset'])
            raw = coordinates.read(row['bytes'])
            n, missing = compare_coordinates(raw, row['model'])
            c = counts[row['source']]
            c['models'] += 1
            c['residues'] += n
            c['missing_backbone_residues'] += missing
            c['source_cif_bytes'] += Path(row['model']['path']).stat().st_size
    if sha(job_path) != job_hash:
        raise ValueError('Geometry job changed during readback')
    result = dict(status='complete_original_cif_to_native_float32_coordinate_shard',
                  job=str(job_path), job_sha256=job_hash, counts=counts)
    with output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    return result


def run_createdb(command, directory):
    directory = Path(directory)
    begin = time.monotonic()
    with (directory / 'createdb.log').open('x') as log:
        native = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                  env=dict(os.environ, CUDA_VISIBLE_DEVICES=''))
        process = psutil.Process(native.pid)
        identity = dict(pid=process.pid, created=process.create_time(), command=command)
        with (directory / 'createdb-original-process.json').open('x') as handle:
            json.dump(identity, handle, indent=2)
            handle.write('\n')
        code = native.wait()
    result = dict(**identity, exit_code=code, wall_seconds=time.monotonic()-begin)
    with (directory / 'createdb-execution.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    if code:
        raise subprocess.CalledProcessError(code, command)
    return result


def geometry_jobs(prefix, model_rows, output, per_shard):
    lookup = {}
    with Path(str(prefix)+'.lookup').open() as handle:
        for line in handle:
            key, name, _ = line.rstrip('\n').split('\t')
            if name in lookup:
                raise ValueError('Repeated native model alias')
            lookup[name] = int(key)
    positions = index(Path(str(prefix)+'_ca.index'))
    if len(lookup) != len(model_rows):
        raise ValueError('Native model cardinality differs')
    output.mkdir()
    jobs, buffer = [], []
    def emit():
        job = output / ('J%07d.jsonl' % len(jobs))
        with job.open('x') as handle:
            for row in buffer:
                handle.write(json.dumps(row, separators=(',', ':'))+'\n')
        jobs.append((str(job), sha(job)))
        buffer.clear()
    for row in model_rows:
        name = Path(row['model']['path']).stem
        offset, size = positions[lookup.pop(name)]
        buffer.append(dict(**row, offset=offset, bytes=size))
        if len(buffer) == per_shard:
            emit()
    if buffer:
        emit()
    if lookup:
        raise ValueError('Foreign native model alias')
    return jobs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    parser.add_argument('--receipt', required=True, type=Path)
    args = parser.parse_args()
    assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text())
    plan_hash = sha(args.plan)
    pins = dict(plan['pins'])
    verify(pins)
    closure = json.loads(Path(plan['atlas_closure']).read_text())
    assert closure['status'] == 'complete_verified_full_prediction_atlas_availability_union'
    assert closure['model_counts'] == plan['expected_models_by_source']
    fixture = json.loads(Path(plan['fixture']).read_text())
    assert fixture['status'] == 'passed_mixed_source_native_database_and_geometry_controls'
    assert fixture['builder_sha256'] == sha(Path(__file__))
    assert fixture['foldseek_sha256'] == sha(plan['foldseek'])
    output = Path(plan['output'])
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    def state(stage, **values):
        temporary = output / 'state.tmp'
        temporary.write_text(json.dumps(dict(stage=stage, **values), indent=2)+'\n')
        temporary.replace(output / 'state.json')
    state('enrolling_full_source_model_grid')
    rows, aliases = [], set()
    counts = {source: dict(models=0, residues=0, source_cif_bytes=0)
              for source in ('AFDB', 'ESMFold')}
    with Path(plan['models']).open() as source, (output / 'model_paths.tsv').open('x') as paths:
        for line in source:
            row = model_record(json.loads(line))
            unique_alias(row, aliases)
            model = row['model']
            path = Path(model['path'])
            if not path.is_file():
                raise FileNotFoundError(path)
            paths.write(str(path.resolve())+'\n')
            c = counts[row['source']]
            c['models'] += 1
            c['residues'] += model['length']
            c['source_cif_bytes'] += path.stat().st_size
            rows.append(row)
            if len(rows) % 100000 == 0:
                state('enrolling_full_source_model_grid', models=len(rows), counts=counts)
    del aliases
    for source in counts:
        assert counts[source]['models'] == plan['expected_models_by_source'][source]
        assert counts[source]['residues'] == plan['expected_residues_by_source'][source]
    assert sha(args.plan) == plan_hash
    verify(pins)
    prefix = output / 'structures'
    command = [plan['foldseek'], 'createdb', str(output / 'model_paths.tsv'), str(prefix),
               '--threads', str(plan['resources']['cpu']), '--gpu', '0',
               '--mask-bfactor-threshold', '70', '--coord-store-mode', '1']
    state('native_createdb_running', counts=counts, command=command)
    native = run_createdb(command, output)
    state('checking_every_native_model_sequence_and_record', counts=counts)
    checked = readback(prefix, [row['model'] for row in rows])
    assert checked == dict(models=sum(c['models'] for c in counts.values()),
                           residues=sum(c['residues'] for c in counts.values()))
    coordinate_path = Path(str(prefix)+'_ca')
    coordinate_hash = sha(coordinate_path)
    jobs = geometry_jobs(prefix, rows, output / 'geometry-jobs', plan['models_per_geometry_shard'])
    del rows
    proofs = output / 'geometry-proofs'
    proofs.mkdir()
    totals = {source: dict(models=0, residues=0, source_cif_bytes=0,
                          missing_backbone_residues=0) for source in counts}
    state('checking_every_original_cif_and_native_coordinate', shards=0, total_shards=len(jobs), counts=totals)
    receipts = []
    with ProcessPoolExecutor(max_workers=plan['resources']['cpu'],
                             mp_context=multiprocessing.get_context('spawn')) as pool:
        futures = {pool.submit(geometry_shard, path, digest, str(coordinate_path),
                    str(proofs / (Path(path).stem+'.json'))): path for path, digest in jobs}
        for future in as_completed(futures):
            result = future.result()
            for source, values in result['counts'].items():
                for key, value in values.items():
                    totals[source][key] += value
            receipts.append(proofs / (Path(result['job']).stem+'.json'))
            state('checking_every_original_cif_and_native_coordinate', shards=len(receipts),
                  total_shards=len(jobs), counts=totals)
    assert sha(coordinate_path) == coordinate_hash
    for source, values in counts.items():
        assert all(totals[source][key] == value for key, value in values.items())
    assert sha(args.plan) == plan_hash
    verify(pins)
    artifacts = {str(path): sha(path) for path in output.iterdir()
                 if path.is_file() and path.name not in ('state.json', 'state.tmp')}
    for path, digest in jobs:
        bind(artifacts, path, digest)
    for path in receipts:
        bind(artifacts, path)
    archive = output / 'source-and-output-bindings.json'
    bindings = dict(pins)
    for path, digest in artifacts.items():
        bind(bindings, path, digest)
    bind(bindings, args.plan)
    archive.write_text(json.dumps(bindings, indent=2)+'\n')
    verify(bindings)
    result = dict(status='complete_full_prediction_atlas_foldseek_database_pending_independent_execution_and_database_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), **checked, counts_by_source=totals,
        plan_sha256=plan_hash, atlas_closure_sha256=sha(plan['atlas_closure']),
        native_createdb_exit_code=native['exit_code'], native_createdb_wall_seconds=native['wall_seconds'],
        geometry_shards=len(jobs), source_and_output_bindings=str(archive),
        source_and_output_bindings_sha256=sha(archive), bound_files=len(bindings),
        source_hashes={path: digest for path, digest in pins.items()},
        coordinate_storage='native float32, every value exactly matches original CIF after float32 conversion',
        seed_mask_threshold=70, seed_mask_is_not_confidence_qualification=True,
        all_sources_and_alternatives_retained=True, original_3di_values_independently_reconstructed=False,
        independent_whole_database_reader_complete=False, scientific_eligibility=False,
        gpu=False, new_predictions=0, new_cost_usd=0,
        scope='All2,961,055source models encoded without dropping either predictor or the unlinked '
              'alternative. Every lookup/AA hash/3Di length and alphabet/coordinate dimension and '
              'finiteness checked; all original CIFs freshly hashed and every native C-alpha '
              'compared exactly after float32 conversion. No cross-source selection, calibrated '
              'confidence/PAE acceptance, independent3Di reconstruction, structural clustering, '
              'homology or evolutionary result. Original execution closure and full independent '
              'database readback remain required.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    state('complete_pending_independent_execution_and_database_readback', counts=totals)
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
