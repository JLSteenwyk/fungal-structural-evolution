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


def main():
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
        pae_recovery=json.loads(recovery.read_text()) if recovery.exists() else dict(state='not_started_or_unavailable'),
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
