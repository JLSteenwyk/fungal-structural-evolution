#!/usr/bin/env python3
"""Measure host and full-atlas storage resources immediately before clustering.

This creates a measurement receipt only.  It neither accepts prerequisite
scientific gates nor launches Foldseek.
"""
from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path


def read_meminfo(path: Path = Path('/proc/meminfo')) -> dict[str, int]:
    values: dict[str, int] = {}
    for line in path.read_text().splitlines():
        fields = line.replace(':', '').split()
        # /proc/meminfo also contains unitless HugePages counters; they are not
        # byte capacities and are deliberately outside this resource receipt.
        if len(fields) != 3 or fields[2] != 'kB':
            continue
        values[fields[0]] = int(fields[1]) * 1024
    if 'MemAvailable' not in values or 'MemTotal' not in values:
        raise ValueError('MemAvailable or MemTotal is absent from meminfo')
    return values


def disk_bytes(path: Path) -> dict[str, int]:
    stat = os.statvfs(path)
    block = stat.f_frsize
    return {
        'total_bytes': stat.f_blocks * block,
        'available_bytes': stat.f_bavail * block,
        'free_bytes': stat.f_bfree * block,
    }


def atlas_artifacts(database: Path) -> tuple[dict[str, int], int]:
    required = ('structures', 'structures_ca', 'structures_h', 'structures_ss', 'structures.lookup')
    result: dict[str, int] = {}
    for name in required:
        path = database / name
        if not path.is_file() or path.stat().st_size <= 0:
            raise ValueError(f'required database artifact absent or empty: {path}')
        result[name] = path.stat().st_size
    return result, sum(result.values())


def lookup_count(path: Path) -> int:
    count = 0
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if len(line.rstrip('\n').split('\t')) != 3:
                raise ValueError(f'invalid lookup row {number}')
            count += 1
    return count


def measure(database: Path, output_parent: Path, expected_models: int, disk_multiplier: float,
            min_available_memory_gib: float) -> dict:
    artifacts, protected_bytes = atlas_artifacts(database)
    models = lookup_count(database / 'structures.lookup')
    if models != expected_models:
        raise ValueError(f'lookup models {models} differs from expected {expected_models}')
    output_parent.mkdir(parents=True, exist_ok=True)
    disk = disk_bytes(output_parent)
    mem = read_meminfo()
    required_free = int(protected_bytes * disk_multiplier)
    required_mem = int(min_available_memory_gib * 2**30)
    passed = disk['available_bytes'] >= required_free and mem['MemAvailable'] >= required_mem
    return {
        'schema_version': 1,
        'status': ('passed_full_atlas_clustering_resource_preflight' if passed
                   else 'failed_full_atlas_clustering_resource_preflight'),
        'measured_utc': datetime.now(timezone.utc).isoformat(),
        'database': str(database.resolve()),
        'output_parent': str(output_parent.resolve()),
        'expected_models': expected_models,
        'observed_lookup_models': models,
        'protected_database_artifacts_bytes': artifacts,
        'protected_database_bytes': protected_bytes,
        'disk': {**disk, 'required_available_bytes': required_free,
                 'required_multiplier_of_protected_database': disk_multiplier},
        'memory': {'total_bytes': mem['MemTotal'], 'available_bytes': mem['MemAvailable'],
                   'required_available_bytes': required_mem},
        'scientific_eligibility': False,
        'scope': ('Fresh host-resource measurement for a future full-atlas native clustering launch. '
                  'It does not verify coordinate/PAE receipts, database identity hashes, Foldseek behavior, '
                  'structural homology, or biological eligibility.'),
    }


def self_test() -> None:
    import tempfile
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        database = root / 'database'
        database.mkdir()
        for name in ('structures', 'structures_ca', 'structures_h', 'structures_ss'):
            (database / name).write_bytes(b'x')
        (database / 'structures.lookup').write_text('1\ta\t0\n2\tb\t1\n')
        result = measure(database, root / 'output', 2, 1.0, 0.0)
        assert result['observed_lookup_models'] == 2
        assert result['protected_database_bytes'] == 16
        (database / 'structures.lookup').write_text('malformed\n')
        try:
            measure(database, root / 'output', 1, 1.0, 0.0)
        except ValueError:
            pass
        else:
            raise AssertionError('malformed lookup accepted')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--database', type=Path)
    parser.add_argument('--output-parent', type=Path)
    parser.add_argument('--expected-models', type=int, default=2_961_055)
    parser.add_argument('--disk-multiplier', type=float, default=3.0)
    parser.add_argument('--minimum-available-memory-gib', type=float, default=64.0)
    parser.add_argument('--receipt', type=Path)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    self_test()
    if args.self_test:
        print(json.dumps({'status': 'passed_full_atlas_clustering_resource_preflight_controls',
                          'checks': ['lookup_count', 'malformed_lookup_rejected']}, indent=2))
        return
    if not args.database or not args.output_parent or not args.receipt:
        parser.error('--database, --output-parent, and --receipt are required')
    if args.receipt.exists():
        raise FileExistsError(f'fresh receipt path required: {args.receipt}')
    if args.disk_multiplier <= 0 or args.minimum_available_memory_gib < 0:
        raise ValueError('resource thresholds must be non-negative and disk multiplier positive')
    result = measure(args.database, args.output_parent, args.expected_models,
                     args.disk_multiplier, args.minimum_available_memory_gib)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, indent=2, sort_keys=True))
    if result['status'].startswith('failed_'):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
