#!/usr/bin/env python3
"""Observe original numerical/source handles and available receipts only."""
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
    pp = Path('metadata/full_weighted_covariance_qualification_plan_20261004_v1.json')
    plan = json.loads(pp.read_text()); verify(plan['pins'])
    inventory = json.loads(Path('metadata/full_weighted_covariance_qualification_launches_20261004_v1.json').read_text())
    assert inventory['source_plan_sha256'] == sha(pp)
    limits = {'cpu.max':'200000 100000','memory.max':str(32 * 2**30),'memory.swap.max':'0'}
    handles = [observe(p, limits) for p in inventory['launches']]
    source_inventory = json.loads(Path('metadata/full_weighted_covariance_source_census_launches_20261004_v1.json').read_text())
    dependency = observe(source_inventory['launches'][2])
    root = Path(plan['output']); receipts = {}
    for name in ['receipt.json','readback.json']:
        p = root / name; record = dict(present=p.exists())
        if p.exists():
            d = json.loads(p.read_text()); record.update(status=d['status'],sha256=sha(p),
                cohorts=d['cohorts'],numerical_audit_rows=d['numerical_audit_rows'],setting_audit_links=d['setting_audit_links'])
        receipts[name] = record
    completed = {}
    for label, p in [('source',plan['source_census_completion']),('numerical',plan['completion'])]:
        p = Path(p); record = dict(path=str(p),present=p.exists())
        if p.exists():
            d = json.loads(p.read_text()); assert sha(d['full_hash_archive']) == d['full_hash_archive_sha256']
            record.update(status=d['status'],sha256=sha(p),bound_source_hashes=d['bound_source_hashes'],exact_process_journals_checked=d['exact_process_journals_checked'])
        completed[label] = record
    launch = json.loads(Path(inventory['launches'][0]).read_text())
    raw = subprocess.check_output(['journalctl','--user','-u',launch['unit'],'-o','json','--no-pager'])
    events = [json.loads(l) for l in raw.splitlines()]
    progress = [r['MESSAGE'] for r in events if r.get('MESSAGE','').startswith('fresh_weighted_covariance_cohorts ')]
    value = dict(status='verified_original_full_four_control_numerical_runtime',checked_utc=datetime.now(timezone.utc).isoformat(),
        source_plan=str(pp),source_plan_sha256=sha(pp),frozen_pins_checked=len(plan['pins']),original_handles=handles,
        original_source_closure_dependency=dependency,expected=plan['expected'],output_root_present=root.exists(),
        cohort_audit_segment_files=len(list((root / 'cohorts').glob('*.audits.jsonl.gz'))),
        cohort_link_segment_files=len(list((root / 'cohorts').glob('*.links.tsv.gz'))),
        receipt_states=receipts,completion_gates=completed,last_progress_message=progress[-1] if progress else None,
        current_failure_files=[str(p) for p in (root / 'failures').rglob('*') if p.is_file()],
        working_model_fits_launched=False,gpu=False,new_cost_usd=0,all_eight_aims_incomplete=True,
        scope='Original identities/journals and actual caps, frozen pins, available receipts/segments and closure archive hashes only. Partial segment files are not completed audits. No independent full output/source replay, native fit, new numerical qualification claim, restart or biological acceptance.')
    with a.output.open('x') as f:json.dump(value,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in value.items() if k not in ['original_handles','original_source_closure_dependency']}),flush=True)


if __name__ == '__main__':main()
