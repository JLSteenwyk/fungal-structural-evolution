#!/usr/bin/env python3
"""Check mixed-source native conversion and original-CIF geometry contracts."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from build_full_atlas_foldseek_database_v1 import (
    compare_coordinates, model_record, run_createdb, unique_alias)
from build_whole_proteome_foldseek_database import index, readback
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--foldseek', required=True)
    parser.add_argument('--models', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--receipt', required=True, type=Path)
    args = parser.parse_args()
    assert not args.output.exists() and not args.receipt.exists()
    pins = {}
    for path in [args.models, args.foldseek, Path(__file__),
                 'scripts/build_full_atlas_foldseek_database_v1.py',
                 'scripts/build_whole_proteome_foldseek_database.py',
                 'scripts/readback_full_atlas_coordinate_profiles_v1.py',
                 'metadata/full_atlas_profile_readback_fixture_20261005_v1.json',
                 'metadata/whole_proteome_foldseek_fixture_checks.json']:
        bind(pins, path)
    selected = {'AFDB': [], 'ESMFold': []}
    with args.models.open() as handle:
        for line in handle:
            raw = json.loads(line)
            if len(selected[raw['source']]) < 2:
                selected[raw['source']].append(model_record(raw))
            if all(len(rows) == 2 for rows in selected.values()):
                break
    assert all(len(rows) == 2 for rows in selected.values())
    rows = selected['AFDB'] + selected['ESMFold']
    aliases = set()
    for row in rows:
        unique_alias(row, aliases)
        bind(pins, row['model']['path'], row['model']['sha256'])
    verify(pins)
    args.output.mkdir(parents=True)
    paths = args.output / 'model_paths.tsv'
    paths.write_text(''.join(str(Path(row['model']['path']).resolve())+'\n' for row in rows))
    prefix = args.output / 'structures'
    command = [args.foldseek, 'createdb', str(paths), str(prefix), '--threads', '1',
               '--gpu', '0', '--mask-bfactor-threshold', '70', '--coord-store-mode', '1']
    native = run_createdb(command, args.output)
    checked = readback(prefix, [row['model'] for row in rows])
    lookup = {}
    with Path(str(prefix)+'.lookup').open() as handle:
        for line in handle:
            key, name, _ = line.rstrip('\n').split('\t')
            lookup[name] = int(key)
    coordinates = index(Path(str(prefix)+'_ca.index'))
    native_buffers = []
    with Path(str(prefix)+'_ca').open('rb') as handle:
        for row in rows:
            offset, size = coordinates[lookup[Path(row['model']['path']).stem]]
            handle.seek(offset)
            raw = handle.read(size)
            assert compare_coordinates(raw, row['model'])[0] == row['model']['length']
            native_buffers.append(raw)
    controls = []
    def rejected(label, operation):
        try:
            operation()
        except (ValueError, RuntimeError):
            controls.append(label)
        else:
            raise AssertionError('Invalid native/source contract accepted: ' + label)
    row = rows[0]
    rejected('duplicate_source_alias', lambda: unique_alias(row, set(aliases)))
    rejected('unknown_prediction_source', lambda: model_record(dict(source='Other', model=row['model'])))
    rejected('boolean_model_length', lambda: model_record(dict(row, model=dict(row['model'], length=True))))
    rejected('invalid_source_digest', lambda: model_record(dict(row, model=dict(row['model'], sha256='bad'))))
    rejected('newline_path_list_injection', lambda: model_record(dict(row, model=dict(row['model'], path='bad\ninput.cif'))))
    bad_models = [dict(rows[0]['model'], sequence_sha256='0'*64)] + [r['model'] for r in rows[1:]]
    rejected('native_sequence_source_disagreement', lambda: readback(prefix, bad_models))
    raw, model = native_buffers[0], rows[0]['model']
    floats = np.frombuffer(raw[:-1], dtype=np.float32).copy()
    floats[0] += 0.5
    rejected('changed_finite_native_coordinate', lambda: compare_coordinates(floats.tobytes()+b'\0', model))
    floats[0] = np.nan
    rejected('nonfinite_native_coordinate', lambda: compare_coordinates(floats.tobytes()+b'\0', model))
    rejected('truncated_native_coordinate_record', lambda: compare_coordinates(raw[:-4], model))
    rejected('foreign_native_record_terminator', lambda: compare_coordinates(raw[:-1]+b'X', model))
    rejected('changed_original_cif_digest', lambda: compare_coordinates(raw, dict(model, sha256='0'*64)))
    for path in args.output.iterdir():
        if path.is_file():
            bind(pins, path)
    verify(pins)
    result = dict(status='passed_mixed_source_native_database_and_geometry_controls',
        checked_utc=datetime.now(timezone.utc).isoformat(), models_by_source={s: len(v) for s, v in selected.items()},
        **checked, native_exit_code=native['exit_code'], native_command=command,
        builder_sha256=sha('scripts/build_full_atlas_foldseek_database_v1.py'),
        foldseek_sha256=sha(args.foldseek), every_fixture_coordinate_exact_float32=True,
        rejection_controls=controls, source_hashes=pins,
        native_outputs_retained=True, scientific_eligibility=False, gpu=False,
        scope='Four native software compatibility controls covering both source representations. '
              'All original CIF atom validation and every native C-alpha/AA/3Di record contract '
              'checked; eleven corruption/identity controls rejected without altering original '
              'source or positive native output. This is not a sampling pilot, full atlas '
              'conversion, independent3Di reconstruction or biological qualification.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
