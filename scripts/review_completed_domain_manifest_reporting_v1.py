#!/usr/bin/env python3
"""Retain a complete manifest whose original wrapper failed while writing its report."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text())
    pins = dict(plan['pins'])
    verify(pins)
    cfg = json.loads(Path(plan['configuration']).read_text())
    original = json.loads(Path(plan['producer_receipt']).read_text())
    payload = json.loads(Path(plan['original_tool']).read_text())
    launch = json.loads(Path(plan['launch']).read_text())
    assert payload['original_tool_session_id'] == payload['initial']['session_id'] == launch['original_tool_session_id']
    assert payload['terminal']['exit_code'] == 1
    assert launch['invocation_id'] == cfg['invocation_id']
    assert launch['pid'] == cfg['wrapper']['pid'] and launch['created'] == cfg['wrapper']['created']
    assert launch['cmdline'] == cfg['wrapper']['cmdline']
    assert cfg['invocation_id'] in payload['initial']['output']
    assert plan['unit'] in payload['initial']['output']
    rows = [json.loads(line) for line in subprocess.check_output(
        ['journalctl', '--user', '-u', plan['unit'], '-o', 'json', '--no-pager'], text=True).splitlines()]
    rows = [r for r in rows if cfg['invocation_id'] in (r.get('_SYSTEMD_INVOCATION_ID'), r.get('USER_INVOCATION_ID'))]
    exact = [r for r in rows if r.get('_PID') == str(cfg['wrapper']['pid'])
             and r.get('_CMDLINE') == ' '.join(cfg['wrapper']['cmdline'])]
    assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=cfg['wrapper'], invocation_id=cfg['invocation_id'])
    messages = [r['MESSAGE'] for r in exact]
    assert sum('Traceback (most recent call last):' in m for m in messages) == 1
    assert any("with a.execution.open('x')" in m for m in messages)
    assert messages[-1] == "FileExistsError: [Errno 17] File exists: '"+plan['colliding_execution_path']+"'"
    assert sum('Main process exited, code=exited, status=1/FAILURE' in r['MESSAGE'] for r in rows) == 1
    assert sum(r.get('CPU_USAGE_NSEC') is not None for r in rows) == 1
    assert Path(plan['colliding_execution_path']).is_dir()
    assert not Path(plan['missing_execution_receipt']).exists()
    assert Path(plan['native_stderr']).read_text() == ''
    assert original['status'] == 'completed_full_refreshed_domain_manifest_pending_independent_readback'
    raw_path = Path(plan['manifest'])/'receipt.json'
    raw = json.loads(raw_path.read_text())
    assert sha(raw_path) == original['raw_receipt_sha256']
    assert raw['status'] == 'complete_all_candidate_domain_extraction_manifest'
    for key in ('source_models', 'unique_intervals', 'boundary_links', 'candidate_model_hit_pairs', 'unique_interval_residues'):
        assert raw[key] == original[key]
    for name, digest in raw['artifacts'].items():
        bind(pins, Path(plan['manifest'])/name, digest)
    for mapping in (cfg['source_hashes'], original['source_hashes'], launch['source_hashes']):
        for path, digest in mapping.items():
            bind(pins, path, digest)
    journal = Path(plan['journal'])
    with journal.open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True)+'\n')
    for path in (args.plan, raw_path, journal, Path(__file__)):
        bind(pins, path)
    verify(pins)
    result = {k: original[k] for k in ('source_models', 'unique_intervals', 'boundary_links', 'candidate_model_hit_pairs', 'unique_interval_residues', 'raw_receipt_sha256')}
    result.update(status='retained_full_domain_manifest_after_wrapper_reporting_failure_pending_independent_readback',
                  checked_utc=datetime.now(timezone.utc).isoformat(), source_hashes=pins,
                  scientific_eligibility=False, gpu=False, new_predictions=0,
                  original_tool_session_id=payload['original_tool_session_id'], original_tool_terminal_exit_code=1,
                  original_native_terminal_exit_code=None, original_native_exit_proof_available=False,
                  native_work_repeated=False, original_wrapper_failure_preserved=True,
                  scope='The original tool wait failed at the final wrapper reporting write because '
                        'its receipt path was also its log directory. Exact original journal, '
                        'identity, empty native stderr and complete bound producer output retained. '
                        'No original native exit code is reconstructed or assumed. A new full '
                        'independent SQL-union readback and its actual original wait must close '
                        'before these saved intervals can enter extraction; no native rerun.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
