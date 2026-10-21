#!/usr/bin/env python3
"""Report live execution state without changing any project artifact."""
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def unit(name):
    result = subprocess.run(
        ['systemctl', '--user', 'show', name, '-p', 'ActiveState', '-p', 'SubState', '-p', 'MainPID'],
        text=True, capture_output=True, check=False)
    fields = dict(line.split('=', 1) for line in result.stdout.splitlines() if '=' in line)
    return dict(unit=name, present=result.returncode == 0 and bool(fields), **fields)


def recovery_with_projection(value):
    """Add a non-binding linear queue projection when state counters permit it."""
    if not isinstance(value, dict):
        return value
    counts = value.get('counts')
    expected = value.get('expected_models')
    elapsed = value.get('elapsed_seconds')
    if not isinstance(counts, dict) or not isinstance(expected, int) or expected < 1:
        return value
    attempted = counts.get('attempted')
    if (not isinstance(attempted, int) or attempted < 1 or attempted > expected
            or not isinstance(elapsed, (int, float)) or elapsed <= 0):
        return value
    rate = attempted / elapsed
    remaining = expected - attempted
    value = dict(value)
    value['linear_queue_projection_advisory'] = {
        'attempts_per_second': rate,
        'remaining_queue_records': remaining,
        'remaining_seconds_at_observed_average': remaining / rate,
        'remaining_days_at_observed_average': remaining / rate / 86400,
        'assumptions': ('Linear extrapolation from all completed recovery attempts so far, including failures. '
                        'It excludes any later serial retry, readback, lineage closure, outages, and changes in '
                        'remote-service behavior; it is not a completion guarantee.'),
    }
    return value


def self_test():
    projected = recovery_with_projection({
        'counts': {'attempted': 50}, 'expected_models': 100, 'elapsed_seconds': 25,
    })
    assert projected['linear_queue_projection_advisory']['attempts_per_second'] == 2
    assert projected['linear_queue_projection_advisory']['remaining_seconds_at_observed_average'] == 25
    assert 'linear_queue_projection_advisory' not in recovery_with_projection({
        'counts': {'attempted': 0}, 'expected_models': 100, 'elapsed_seconds': 25,
    })


def main():
    self_test()
    recovery = ROOT / 'results/structures/full-atlas-missing-pae-recovery-20261021-v1/state.json'
    pmsf = ROOT / 'results/phylogeny/native-taxon-PMSF-refits-20261001-v1'
    pmsf_completed = len(list(pmsf.glob('*/receipt.json')))
    repeat_audits = sorted((ROOT / 'metadata').glob('repeatmodeler_full_panel_progress_audit_*.json'))
    repeat = None
    if repeat_audits:
        payload = json.loads(repeat_audits[-1].read_text())
        repeat = dict(audit=str(repeat_audits[-1].relative_to(ROOT)), counts=payload.get('counts', {}))
    scheduler = subprocess.run(['scontrol', 'ping'], text=True, capture_output=True, check=False)
    value = dict(
        status='live_pipeline_status_advisory',
        pae_recovery=(recovery_with_projection(json.loads(recovery.read_text()))
                      if recovery.exists() else dict(state='not_started_or_unavailable')),
        pmsf=dict(completed_receipts=pmsf_completed, expected_fits=16),
        repeatmodeler=repeat,
        scheduler=dict(returncode=scheduler.returncode, message=(scheduler.stdout + scheduler.stderr).strip()),
        services=[unit(name) for name in (
            'fungal-full-atlas-missing-pae-recovery-20261021-v1.service',
            'fungal-full-atlas-missing-pae-recovery-readback-20261021-v1.service',
            'fungal-full-atlas-pae-failure-retry-20261021-v1.service',
            'fungal-full-atlas-pae-failure-retry-readback-20261021-v1.service',
            'fungal-full-atlas-pae-recovery-lineage-20261021-v2.service',
            'fungal-repeatmodeler-recovery-submit-20261021-v1.service',
            'fungal-repeatmodeler-recovery-audit-20261021-v1.service',
            'fungal-species-taxon-pmsf-resume-20261021-v1.service',
            'fungal-species-taxon-pmsf-resume-collection-20261021-v1.service',
            'fungal-species-taxon-pmsf-resume-collection-closure-20261021-v1.service')],
        scope='Read-only runtime observation. It establishes neither completed scientific stages nor biological inference.')
    print(json.dumps(value, indent=2))


if __name__ == '__main__':
    main()
