#!/usr/bin/env python3
"""Observe original design/qualification/timing/publication handles without restarting them."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal, live_record
from run_ortholog_pair_guide_comparison import sha


PLAN_PATHS = [
    'metadata/full_expanded_model_designs_plan_20261002_v2.json',
    'metadata/full_uniform_covariance_qualification_plan_20261002.json',
    'metadata/full_shared_entity_fit_draft_plan_20261002_v2.json',
    'metadata/full_shared_entity_timing_plan_20261002.json',
    'metadata/baliphy_recovery_full_diagnostics_plan_20261002.json',
    'metadata/baliphy_full_recovered_horizon_diagnostic_publication_plan_20261002.json',
]
LAUNCH_STEMS = [
    'full_expanded_model_designs_v2', 'full_expanded_model_designs_v2_readback',
    'full_expanded_model_designs_v2_closure', 'full_uniform_covariance_qualification',
    'full_uniform_covariance_qualification_readback', 'full_uniform_covariance_qualification_closure',
    'full_shared_entity_timing', 'full_shared_entity_timing_readback', 'full_shared_entity_timing_closure',
    'baliphy_recovery_full_diagnostics', 'baliphy_recovery_full_diagnostics_readback',
    'baliphy_recovery_full_diagnostics_closure', 'baliphy_recovered_diagnostics_publication',
    'baliphy_recovered_diagnostics_publication_closure',
]
COMPLETIONS = {
    'design': ('metadata/full_expanded_model_designs_v2_completed_20261002.json', 'complete_verified_full_expanded_model_designs'),
    'qualification': ('metadata/full_uniform_covariance_qualification_completed_20261002.json', 'complete_verified_full_uniform_covariance_qualification'),
    'timing': ('metadata/full_shared_entity_timing_completed_20261002.json', 'complete_verified_full_shared_entity_timing_accounting'),
    'recovery_diagnostics': ('metadata/baliphy_recovery_full_diagnostics_completed_20261002.json', 'complete_verified_full_baliphy_recovery_diagnostics'),
    'recovered_publication': ('metadata/baliphy_full_recovered_horizon_diagnostic_publication_completed_20261002.json', 'complete_verified_full_recovered_horizon_diagnostic_publication'),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    bindings = {str(Path(__file__)): sha(Path(__file__))}
    for path in PLAN_PATHS:
        plan = json.loads(Path(path).read_text())
        bindings[path] = sha(path)
        for p, d in plan['pins'].items():
            assert p not in bindings or bindings[p] == d, p
            bindings[p] = d
    assert all(sha(p) == d for p, d in bindings.items())
    handles = []
    for stem in LAUNCH_STEMS:
        path = 'metadata/' + stem + '_launch_20261002.json'
        record = json.loads(Path(path).read_text())
        record['launch'] = path
        assert sha(record['plan']) == record['plan_sha256']
        process = fingerprint(record)
        if process is None:
            observation = journal_terminal(record)
        else:
            observation = live_record(record, process)
            cgroup = Path('/proc') / str(record['pid']) / 'cgroup'
            try:
                lines = cgroup.read_text().splitlines()
                unified = [line[3:] for line in lines if line.startswith('0::')]
                assert len(unified) == 1
                root = Path('/sys/fs/cgroup') / unified[0].lstrip('/')
                limits = {k: (root / k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
                assert limits == record['actual_cgroup_limits']
                assert limits == {'cpu.max': '200000 100000', 'memory.max': str(32 * 2**30), 'memory.swap.max': '0'}
                observation['current_cgroup_limits'] = limits
            except FileNotFoundError:
                assert fingerprint(record) is None
                observation = journal_terminal(record)
        handles.append(observation)
    stages = {}
    for label, (path, expected) in COMPLETIONS.items():
        if not Path(path).exists():
            stages[label] = dict(path=path, accounting_closure_present=False)
            continue
        record = json.loads(Path(path).read_text())
        assert record['status'] == expected, (path, record['status'])
        assert sha(record['full_hash_archive']) == record['full_hash_archive_sha256']
        stages[label] = dict(path=path, sha256=sha(path), accounting_closure_present=True,
            status=record['status'], full_hash_archive=record['full_hash_archive'],
            full_hash_archive_sha256=record['full_hash_archive_sha256'],
            reported_full_source_bindings=record['bound_source_hashes'],
            reported_original_process_journals=record['exact_process_journals_checked'])
    fit = json.loads(Path(PLAN_PATHS[2]).read_text())
    assert fit['launch_state'] == 'not_launched_or_queued' and not Path(fit['output']).exists()
    result = dict(status='verified_original_full_timing_queue_and_recovered_publication_checkpoint',
        checked_utc=datetime.now(timezone.utc).isoformat(), source_hashes=bindings,
        original_handles=handles, stage_observations=stages,
        live_original_handles=sum('current_cgroup_limits' in r for r in handles),
        original_terminal_success_journals=sum(r['status'].startswith('verified_original_terminal') for r in handles),
        production_fitting_launched=False, gpu_inference_paused=True, all_eight_aims_incomplete=True,
        scope='Small pinned plan/backend sources and compact archive hashes rechecked; every captured original PID/create/command, current live cgroup or actual invocation-linked terminal resource journal observed. Full multi-million-binding scientific archives are not rehashed by this observer, and a missing compact closure is not a process failure. No queue/source mutation, restart, GPU prediction or inferential acceptance.')
    with args.output.open('x') as f:
        f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'original_handles']}))


if __name__ == '__main__':
    main()
