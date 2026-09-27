"""Archive terminal states, pinned provenance and completed accessibility stages.

Run from the project root. This verifies completion artifacts; the separate
full readbacks establish row-level checks, not biological exposure accuracy.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import subprocess


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(output)
    checked = {}

    def check(path, expected):
        actual = sha(path)
        if actual != expected:
            raise ValueError(f'Hash mismatch: {path}')
        checked[str(path)] = actual

    stages = []
    units = {}
    for name in ('accessibility', 'accessibility_projection', 'site_exposure'):
        stem = f'metadata/recovered_afdb_{name}'
        plan_path = Path(stem + '_plan_20260927.json')
        launch_path = Path(stem + '_launch_20260927.json')
        plan = json.loads(plan_path.read_text())
        launch = json.loads(launch_path.read_text())
        check(plan_path, launch['plan_sha256'])
        checked[str(launch_path)] = sha(launch_path)
        raw = subprocess.check_output([
            'systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState',
            '-p', 'SubState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True)
        state = dict(line.split('=', 1) for line in raw.splitlines() if '=' in line)
        if state != {'ActiveState': 'inactive', 'SubState': 'dead',
                     'Result': 'success', 'ExecMainStatus': '0'}:
            raise ValueError(f'Producer not successfully terminal: {state}')
        units[launch['unit']] = state
        for path, expected in plan['pinned_files'].items():
            check(path, expected)
        controller_path = Path(plan['output']) / 'receipt.json'
        controller = json.loads(controller_path.read_text())
        checked[str(controller_path)] = sha(controller_path)
        if controller['config_sha256'] != sha(plan_path):
            raise ValueError('Controller plan mismatch')
        for path, expected in controller['prerequisite_receipts'].items():
            check(path, expected)
        recorded = {s['stage']: s for s in controller['stages']}
        if set(recorded) != {s['name'] for s in plan['stages']}:
            raise ValueError('Stage inventory mismatch')
        for stage in plan['stages']:
            path = Path(stage['receipt'])
            source = recorded[stage['name']]
            if source['receipt'] != str(path):
                raise ValueError('Stage receipt path mismatch')
            check(path, source['receipt_sha256'])
            receipt = json.loads(path.read_text())
            if receipt['status'] != stage['status']:
                raise ValueError(f'Incomplete stage: {path}')
            for key, value in stage.get('expected_counts', {}).items():
                if receipt[key] != value:
                    raise ValueError(f'Count mismatch: {path}: {key}')
            for artifact, expected in receipt.get('artifacts', {}).items():
                check(path.parent / artifact, expected)
            stages.append({'receipt': str(path), 'status': receipt['status'],
                           'expected_counts_verified': stage.get('expected_counts', {})})
    result = {
        'status': 'passed_recovered_accessibility_completion_provenance_check',
        'checked_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'script_sha256': sha(__file__), 'terminal_units': units,
        'stages': stages, 'checked_sha256': checked,
        'scope': 'Terminal success, launch-bound plans, all plan pins, controller prerequisites, stage receipts/counts and declared stage artifacts. Row-level numeric checks are supplied by the separate completed full readbacks. ASA geometry was not independently recalculated; no ancestral exposure, rate effect, or biological core/surface validation is claimed.'}
    with output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(f'Checked {len(stages)} stages and {len(checked)} unique files.')


if __name__ == '__main__':
    main()
