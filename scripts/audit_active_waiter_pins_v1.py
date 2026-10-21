#!/usr/bin/env python3
"""Audit that hash-pinned inputs of queued dependency plans remain unchanged."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile
from typing import Any

DEFAULT_GLOB = '*wait_plan_20261021_v*.json'


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def audit_plan(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f'plan must be a JSON object: {path}')
    pins = value.get('pins')
    if not isinstance(pins, dict):
        raise ValueError(f'plan pins must be an object: {path}')
    missing: list[str] = []
    mismatches: list[str] = []
    malformed: list[str] = []
    for raw_path, expected in sorted(pins.items()):
        if not isinstance(raw_path, str) or not isinstance(expected, str) or len(expected) != 64:
            malformed.append(str(raw_path))
            continue
        target = Path(raw_path)
        if not target.is_file():
            missing.append(raw_path)
        elif sha256(target) != expected:
            mismatches.append(raw_path)
    return {
        'plan': str(path), 'pins': len(pins), 'missing': missing,
        'mismatches': mismatches, 'malformed_pin_records': malformed,
        'passed': not missing and not mismatches and not malformed,
    }


def audit(paths: list[Path]) -> dict[str, Any]:
    if not paths:
        raise ValueError('at least one plan is required')
    plans = [audit_plan(path) for path in paths]
    return {
        'status': ('passed_active_waiter_pin_audit' if all(item['passed'] for item in plans)
                   else 'failed_active_waiter_pin_audit'),
        'plans': plans,
        'scientific_eligibility': False,
        'scope': ('Integrity audit for queued dependency plans only. It verifies recorded input hashes, '
                  'not process liveness, completed analysis receipts, or scientific eligibility.'),
    }


def self_test() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        input_path = root / 'input.txt'
        input_path.write_text('original\n')
        plan = root / 'plan.json'
        plan.write_text(json.dumps({'pins': {str(input_path): sha256(input_path)}}))
        assert audit([plan])['status'] == 'passed_active_waiter_pin_audit'
        input_path.write_text('changed\n')
        result = audit([plan])
        assert result['status'] == 'failed_active_waiter_pin_audit'
        assert result['plans'][0]['mismatches'] == [str(input_path)]
        plan.write_text(json.dumps({'pins': {str(root / 'missing.txt'): '0' * 64}}))
        result = audit([plan])
        assert result['plans'][0]['missing'] == [str(root / 'missing.txt')]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', action='append', type=Path, default=[],
                        help='queued dependency plan to audit; may be repeated')
    parser.add_argument('--metadata-dir', type=Path, default=Path('metadata'))
    parser.add_argument('--glob', default=DEFAULT_GLOB,
                        help='plan glob used when --plan is omitted')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    self_test()
    if args.self_test:
        print(json.dumps({'status': 'passed_active_waiter_pin_audit_controls',
                          'checks': ['matching_pin', 'changed_input_rejected', 'missing_input_rejected']}, indent=2))
        return
    paths = args.plan or sorted(args.metadata_dir.glob(args.glob))
    result = audit(paths)
    print(json.dumps(result, indent=2, sort_keys=True))
    if result['status'].startswith('failed_'):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
