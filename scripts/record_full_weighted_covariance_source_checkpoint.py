#!/usr/bin/env python3
"""Observe original full four-control source jobs without restarting or replaying."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe
from reference_measurement_union_sources import verify


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists()
    pp = Path('metadata/full_weighted_covariance_source_census_plan_20261004_v1.json')
    plan = json.loads(pp.read_text()); verify(plan['pins'])
    ip = Path('metadata/full_weighted_covariance_source_census_launches_20261004_v1.json')
    inventory = json.loads(ip.read_text()); assert inventory['source_plan_sha256'] == sha(pp)
    expected = {'cpu.max': '200000 100000', 'memory.max': str(16 * 2**30), 'memory.swap.max': '0'}
    handles = [observe(p, expected) for p in inventory['launches']]
    launch = json.loads(Path(inventory['launches'][0]).read_text())
    raw = subprocess.check_output(['journalctl', '--user', '-u', launch['unit'], '-o', 'json', '--no-pager'])
    events = [json.loads(l) for l in raw.splitlines()]
    progress = [r for r in events if r.get('MESSAGE', '').startswith('reconstructed_original_weighted_source_cohorts ')]
    root = Path(plan['output']); counts = {}
    for name in ['receipt.json', 'readback.json']:
        p = root / name; entry = dict(present=p.exists())
        if p.exists():
            value = json.loads(p.read_text()); entry.update(sha256=sha(p), status=value['status'],
                cohorts=value['cohorts'], designs=value['designs'], fit_inputs=value['fit_inputs'], settings=value['settings'])
        counts[name] = entry
    cp = Path(plan['completion']); closure = dict(path=str(cp), present=cp.exists())
    if cp.exists():
        c = json.loads(cp.read_text()); assert sha(c['full_hash_archive']) == c['full_hash_archive_sha256']
        closure.update(status=c['status'], sha256=sha(cp), bound_source_hashes=c['bound_source_hashes'],
            exact_process_journals_checked=c['exact_process_journals_checked'])
    value = dict(status='verified_original_full_four_control_source_census_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), source_plan=str(pp), source_plan_sha256=sha(pp),
        frozen_pins_checked=len(plan['pins']), original_handles=handles, expected=plan['expected'],
        output_root_present=root.exists(), receipt_states=counts, completion=closure,
        last_cohort_progress_message=progress[-1]['MESSAGE'] if progress else None,
        journal_snapshot_sha256=__import__('hashlib').sha256(raw).hexdigest(),
        numerical_audits_computed=0, working_model_fits_computed=0, gpu=False, settings_changed=False,
        scope='Original PID/create/command or exact terminal journal, live caps/resources, frozen pins and available receipts. No fresh entire consumed source replay, covariance numerical qualification, restart or posterior/biological acceptance. Running source checks may precede the output root.')
    with a.output.open('x') as f: json.dump(value, f, indent=2); f.write('\n')
    print(json.dumps({k:v for k,v in value.items() if k != 'original_handles'}), flush=True)


if __name__ == '__main__': main()
