#!/usr/bin/env python3
"""Check actual full-extraction custody plus explicit identity/terminal corruptions."""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_completed_domain_coordinate_native_closure_v1 import NAMES, inputs, validate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', required=True, type=Path)
    args = parser.parse_args()
    assert not args.receipt.exists()
    docs = inputs()
    result = validate(docs)
    assert result['original_tool_terminal_exit_code'] is None and result['original_native_exit_code'] == 0
    cases = {
        'native_exit_changed': lambda d: d['execution'].update(exit_code=1),
        'timeout_changed': lambda d: d['execution'].update(timed_out=True),
        'wrapper_identity_changed': lambda d: d['launch'].update(pid=d['launch']['pid'] + 1),
        'invocation_changed': lambda d: d['launch'].update(invocation_id='0' * 32),
        'child_command_changed': lambda d: d['process']['command'].append('--invented'),
        'api_zero_fabricated': lambda d: d['initial'].update(terminal={'exit_code': 0}),
        'api_unavailable_mislabelled': lambda d: d['unavailable'].update(tool_call_failed=False),
        'missing_manager_completion': lambda d: d.update(journal=[r for r in d['journal'] if not r.get('CPU_USAGE_NSEC')]),
        'changed_wrapper_terminal': lambda d: d['journal'][-2].update(MESSAGE='{}'),
        'partial_shard_scope': lambda d: d['raw'].update(shards=976),
        'changed_global_counts': lambda d: d['raw']['counts'].update(intervals=2454564),
    }
    rejected = []
    for name, mutate in cases.items():
        altered = copy.deepcopy(docs)
        mutate(altered)
        try:
            validate(altered)
        except (AssertionError, KeyError, ValueError, TypeError):
            rejected.append(name)
        else:
            raise AssertionError('Altered original custody accepted: ' + name)
    paths = [*map(Path, NAMES.values()), Path(__file__), Path('scripts/record_completed_domain_coordinate_native_closure_v1.py'),
        Path('scripts/ancestral_chain_attempt.py'), Path('scripts/reference_measurement_union_sources.py')]
    ep = Path(NAMES['execution']).with_suffix('')
    paths.extend([ep / 'configuration.json', ep / 'process.json', Path(docs['plan']['output']) / 'receipt.json'])
    receipt = dict(status='passed_domain_native_closure_actual_records_and_identity_rejections',
        checked_utc=datetime.now(timezone.utc).isoformat(), actual_original_native_records_checked=True,
        original_tool_terminal_exit_code=None, original_native_exit_code=0, original_api_unknown_preserved=True,
        full_models=976357, full_intervals=2454565, full_shards=977, rejection_controls=rejected,
        closure_script_sha256=sha('scripts/record_completed_domain_coordinate_native_closure_v1.py'),
        source_hashes={str(p): sha(p) for p in paths}, native_work_repeated=False, scientific_eligibility=False,
        scope='Actual full original producer custody and eleven mutated identity/scope/terminal controls. '
              'API terminal remains unknown. This checks native closure software, not exported atoms.')
    with args.receipt.open('x') as handle:
        json.dump(receipt, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
