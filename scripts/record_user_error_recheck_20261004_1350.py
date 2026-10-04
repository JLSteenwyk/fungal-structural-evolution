#!/usr/bin/env python3
"""Record fresh formatter probes and independently reread existing error evidence."""
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from read_baliphy_scalar_json_v6b import compare_tsv


def main():
    output = Path('metadata/current_analysis_error_recheck_20261004_user_1350.json')
    assert not output.exists()
    pins = {}

    def read(path):
        path = Path(path)
        pins[str(path)] = sha(path)
        return json.loads(path.read_text())

    def verify(bindings):
        for path, digest in bindings.items():
            assert sha(path) == digest, path
            pins[path] = digest

    prefix = 'metadata/baliphy_logging_error_recheck_'
    receipt = read(prefix + '20261004_user_1350.json')
    execution = read(prefix + 'execution_20261004_user_1350.json')
    tool = read('metadata/current_analysis_error_original_tool_payload_20261004_user_1350.json')
    assert tool['exit_code'] == 0 and 'Finished with result: success' in tool['output']
    assert execution['exit_code'] == 0 and not execution['timed_out']
    assert execution['receipt_sha256'] == sha(execution['receipt'])
    assert execution['invocation_id'] in tool['output']
    verify(receipt['source_hashes'])
    verify(execution['source_hashes'])
    verify(execution['artifacts'])
    journal = Path(prefix + 'original_journal_20261004_user_1350.jsonl')
    messages = [json.loads(line)['MESSAGE'] for line in journal.read_text().splitlines()]
    initial, terminal = map(json.loads, messages)
    assert initial['original_wrapper'] == execution['wrapper']
    assert initial['invocation_id'] == execution['invocation_id']
    excluded = {'source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope'}
    assert terminal == {k: v for k, v in execution.items() if k not in excluded}
    pins[str(journal)] = sha(journal)

    qualified = read('metadata/baliphy_scalar_json_v6_software_validation_20261004_v9.json')
    verify(qualified['source_hashes'])
    comparisons = [compare_tsv(item['new_directory']) for item in qualified['paired_prior_checks']]
    assert sum(item['mapped_values_compared'] for item in comparisons) == 2709

    previous = read('metadata/current_analysis_error_recheck_20261004_user_1224.json')
    failures = {p: h for p, h in previous['source_hashes'].items()
                if '/full-baliphy-scalar-v6-short-sampler-' in p and p.endswith('/receipt.json')}
    assert len(failures) == 24
    verify(failures)
    assert all(read(path)['exit_code'] == -11 for path in failures)

    stack = read('metadata/scalar_native_stack_comparison_terminal_review_20261004_v3.json')
    verify(stack['source_hashes'])
    assert stack['native_signal_from_original_debugger_text'] == 'SIGKILL'
    assert stack['last_saved_iteration'] == 9 and not stack['planned_horizon_completed']
    assert not stack['standalone_stack_change_validated_as_fix']
    pins[str(Path(__file__))] = sha(__file__)
    verify(pins.copy())
    result = dict(
        status='completed_fresh_user_error_recheck',
        checked_utc=datetime.now(timezone.utc).isoformat(),
        installed_formatter_error_reproduced=receipt['numeric_probes']['original']['error_reproduced'],
        incorrect_native_constants=len(receipt['numeric_probes']['original']['incorrect_constant_indices']),
        example=dict(expected=2.34e-10, observed=receipt['numeric_probes']['original']['values'][0]),
        corrected_native_roundtrips=receipt['numeric_probes']['cjson']['native_exact_roundtrips'],
        corrected_v6_saved_mapped_values_freshly_compared=2709,
        original_native_failures_freshly_hash_verified=24,
        original_failures_family='OG0000972',
        original_failure_signal='SIGSEGV',
        stack_comparison_native_signal='SIGKILL',
        stack_comparison_saved_last_iteration=9,
        stack_comparison_planned_iterations=20,
        validated_sampler_fix=False,
        fresh_formatter_execution_invocation_id=execution['invocation_id'],
        fresh_formatter_original_tool_exit_code=tool['exit_code'],
        whole_wrapper_journal_messages_matched=True,
        source_hashes=pins,
        production_restarted=False, installed_software_changed=False,
        new_mcmc_runs=0, scientific_eligibility=False,
        scope='Two fresh bounded native formatter probes reproduce the original numeric bug '
              'and confirm the corrected twelve-value path. Existing corrected scalar fixtures '
              'are independently reread. Original failed sampler receipts and closed diagnostic '
              'artifacts are freshly hash-verified; no sampler crash is rerun or repaired.')
    with output.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
