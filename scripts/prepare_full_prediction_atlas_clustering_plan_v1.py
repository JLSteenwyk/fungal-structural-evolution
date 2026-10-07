#!/usr/bin/env python3
"""Bind validated full-atlas inputs into one immutable native-clustering plan."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile


EXPECTED_MODELS = 2_961_055


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            value.update(block)
    return value.hexdigest()


def load_object(path: Path) -> dict:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f'JSON object required: {path}')
    return value


def completed_receipt(path: Path) -> tuple[dict, str]:
    value = load_object(path)
    status = value.get('status')
    if not isinstance(status, str) or not (status.startswith('completed_') or status.startswith('passed_')):
        raise ValueError(f'completed or passed receipt required: {path}')
    return value, status


def lookup_count(path: Path) -> int:
    identifiers: set[str] = set()
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 3 or not fields[1] or fields[1] in identifiers:
                raise ValueError(f'invalid or repeated lookup identity at row {number}')
            identifiers.add(fields[1])
    return len(identifiers)


def database_files(database: Path) -> dict[str, str]:
    required = ('structures', 'structures_ca', 'structures_h', 'structures_ss', 'structures.lookup')
    result: dict[str, str] = {}
    for name in required:
        path = database / name
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f'required database file absent or empty: {path}')
    for path in sorted(database.iterdir()):
        if path.is_file() and not path.name.endswith('.log') and path.name not in {'state.json'}:
            result[path.name] = sha256(path)
    return result


def build_plan(args: argparse.Namespace) -> dict:
    database = args.database.resolve()
    lookup = database / 'structures.lookup'
    models = lookup_count(lookup)
    if models != args.expected_models:
        raise ValueError(f'database lookup has {models}, expected {args.expected_models}')
    artifacts = database_files(database)
    source_map = args.source_identity_map.resolve()
    identity_receipt = load_object(args.identity_receipt)
    if identity_receipt.get('status') != 'completed_full_atlas_cluster_identity_map':
        raise ValueError('completed identity-map receipt required')
    if Path(identity_receipt.get('output', '')).resolve() != source_map:
        raise ValueError('identity-map receipt output differs from source map')
    if identity_receipt.get('output_sha256') != sha256(source_map):
        raise ValueError('identity-map receipt hash differs from source map')
    if identity_receipt.get('models') != models:
        raise ValueError('identity-map receipt model count differs from lookup')
    if not source_map.is_file() or source_map.stat().st_size == 0:
        raise ValueError(f'source-identity map absent or empty: {source_map}')
    producer, producer_status = completed_receipt(args.database_receipt)
    if producer.get('models') != models:
        raise ValueError('database producer model count differs from lookup')
    reader, _ = completed_receipt(args.database_readback)
    if reader.get('producer_receipt_sha256') != sha256(args.database_receipt):
        raise ValueError('database independent readback is not bound to database producer receipt')
    resource = load_object(args.resource_preflight)
    if resource.get('status') != 'passed_full_atlas_clustering_resource_preflight':
        raise ValueError('passed fresh clustering resource preflight required')
    if resource.get('observed_lookup_models') != models:
        raise ValueError('resource preflight lookup count differs from database')
    observed_database = resource.get('database')
    if observed_database and Path(observed_database).resolve() != database:
        raise ValueError('resource preflight measured another database')
    required_receipts = {}
    pins = {
        str(args.database_receipt): sha256(args.database_receipt),
        str(args.database_readback): sha256(args.database_readback),
        str(args.resource_preflight): sha256(args.resource_preflight),
        str(args.identity_receipt): sha256(args.identity_receipt),
        str(source_map): sha256(source_map),
    }
    for name, path in (
        ('coordinate_profiles', args.coordinate_receipt),
        ('missing_pae', args.pae_receipt),
        ('missing_pae_readback', args.pae_readback_receipt),
    ):
        _, status = completed_receipt(path)
        required_receipts[name] = {'path': str(path), 'sha256': sha256(path), 'status': status}
        pins[str(path)] = sha256(path)
    foldseek = args.foldseek.resolve()
    if not foldseek.is_file():
        raise ValueError(f'Foldseek executable absent: {foldseek}')
    if not args.cluster_argument:
        raise ValueError('one or more explicit --cluster-argument values are required')
    if args.cpu < 1 or args.minimum_free_disk_gib < 1 or args.emergency_free_disk_gib < 1:
        raise ValueError('resource values must be positive')
    if args.emergency_free_disk_gib >= args.minimum_free_disk_gib:
        raise ValueError('emergency disk reserve must be below minimum free disk')
    return {
        'schema_version': 1,
        'status': 'prepared_full_atlas_native_clustering_plan',
        'database': str(database),
        'database_prefix': 'structures',
        'lookup_name': 'structures.lookup',
        'expected_models': models,
        'database_artifacts': artifacts,
        'source_identity_map': str(source_map),
        'source_identity_map_sha256': sha256(source_map),
        'source_identity_receipt': {'path': str(args.identity_receipt), 'sha256': sha256(args.identity_receipt)},
        'database_producer_receipt': {'path': str(args.database_receipt), 'sha256': sha256(args.database_receipt),
                                      'status': producer_status},
        'database_readback_receipt': {'path': str(args.database_readback), 'sha256': sha256(args.database_readback)},
        'required_receipts': required_receipts,
        'resource_preflight': {'path': str(args.resource_preflight), 'sha256': sha256(args.resource_preflight)},
        'pins': pins,
        'foldseek': str(foldseek),
        'cluster_arguments': args.cluster_argument,
        'resources': {'cpu': args.cpu, 'minimum_free_disk_gib': args.minimum_free_disk_gib,
                      'emergency_free_disk_gib': args.emergency_free_disk_gib},
        'output': str(args.output),
        'scientific_eligibility': False,
        'scope': ('Immutable launch plan for a candidate full-atlas structural partition. It does not establish '
                  'structural homology, orthology, confidence calibration, or an evolutionary result.'),
    }


def self_test() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary); db = root / 'db'; db.mkdir()
        for name in ('structures', 'structures_ca', 'structures_h', 'structures_ss'):
            (db / name).write_text('x')
        (db / 'structures.lookup').write_text('1\ta\t0\n2\tb\t1\n')
        source = db / 'source-map.json'; source.write_text('{}\n')
        def receipt(name: str, status: str, **extra) -> Path:
            path = root / name; path.write_text(json.dumps({'status': status, **extra})); return path
        producer = receipt('producer.json', 'completed_database', models=2)
        database_reader = receipt('database-reader.json', 'passed_database_reader', producer_receipt_sha256=sha256(producer))
        coordinate = receipt('coordinate.json', 'passed_coordinate_reader')
        pae = receipt('pae.json', 'completed_pae')
        pae_reader = receipt('pae-reader.json', 'passed_pae_reader')
        resource = receipt('resource.json', 'passed_full_atlas_clustering_resource_preflight',
                           observed_lookup_models=2, database=str(db.resolve()))
        identity = receipt('identity.json', 'completed_full_atlas_cluster_identity_map', output=str(source), output_sha256=sha256(source), models=2)
        args = argparse.Namespace(database=db, database_receipt=producer, database_readback=database_reader,
            source_identity_map=source, identity_receipt=identity, resource_preflight=resource, coordinate_receipt=coordinate,
            pae_receipt=pae, pae_readback_receipt=pae_reader, expected_models=2, foldseek=Path('/bin/true'),
            cluster_argument=['--example'], cpu=1, minimum_free_disk_gib=2, emergency_free_disk_gib=1,
            output=root / 'out')
        assert build_plan(args)['expected_models'] == 2
        pae.write_text(json.dumps({'status': 'running'}))
        try:
            build_plan(args)
        except ValueError:
            pass
        else:
            raise AssertionError('running prerequisite receipt accepted')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--database', type=Path)
    parser.add_argument('--database-receipt', type=Path)
    parser.add_argument('--database-readback', type=Path)
    parser.add_argument('--source-identity-map', type=Path)
    parser.add_argument('--identity-receipt', type=Path)
    parser.add_argument('--resource-preflight', type=Path)
    parser.add_argument('--coordinate-receipt', type=Path)
    parser.add_argument('--pae-receipt', type=Path)
    parser.add_argument('--pae-readback-receipt', type=Path)
    parser.add_argument('--foldseek', type=Path)
    parser.add_argument('--cluster-argument', action='append', default=[])
    parser.add_argument('--cpu', type=int, default=8)
    parser.add_argument('--minimum-free-disk-gib', type=int, default=128)
    parser.add_argument('--emergency-free-disk-gib', type=int, default=64)
    parser.add_argument('--expected-models', type=int, default=EXPECTED_MODELS)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    self_test()
    if args.self_test:
        print(json.dumps({'status': 'passed_full_atlas_clustering_plan_controls',
                          'checks': ['completed_receipt_binding', 'running_receipt_rejected',
                                     'lookup_identity_count', 'resource_preflight_binding']}, indent=2))
        return
    needed = ('database', 'database_receipt', 'database_readback', 'source_identity_map', 'identity_receipt', 'resource_preflight',
              'coordinate_receipt', 'pae_receipt', 'pae_readback_receipt', 'foldseek', 'output', 'plan')
    if any(getattr(args, item) is None for item in needed):
        parser.error('all input paths, --output and fresh --plan are required')
    if args.plan.exists() or args.output.exists():
        raise FileExistsError('fresh plan and output paths are required')
    plan = build_plan(args)
    args.plan.parent.mkdir(parents=True, exist_ok=True)
    args.plan.write_text(json.dumps(plan, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': plan['status'], 'plan': str(args.plan),
                      'models': plan['expected_models'], 'pinned_inputs': len(plan['pins'])}, indent=2))


if __name__ == '__main__':
    main()
