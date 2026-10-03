#!/usr/bin/env python3
"""Observe original full retained timing handles and available closure evidence."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe
from reference_measurement_union_sources import verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    pp = Path('metadata/full_retained_shared_entity_timing_plan_20261003_v1.json')
    plan = json.loads(pp.read_text()); verify(plan['pins'])
    inventory = json.loads(Path(plan['launch_inventory']).read_text())
    assert inventory['source_plan_sha256'] == sha(pp)
    limits = {'cpu.max': '200000 100000', 'memory.max': str(32 * 2**30), 'memory.swap.max': '0'}
    handles = [observe(p, limits) for p in inventory['launches']]
    handles.append(observe(plan['retained_closure_launch']))
    gates = {}
    fit = json.loads(Path(plan['fit_plan']).read_text())
    assert sha(plan['fit_plan']) == plan['fit_plan_sha256']
    for path in [fit['qualification_completion'], fit['retained_completion'], plan['completion']]:
        fp = Path(path); entry = dict(present=fp.exists())
        if fp.exists():
            d = json.loads(fp.read_text())
            assert sha(d['full_hash_archive']) == d['full_hash_archive_sha256']
            entry.update(status=d['status'], sha256=sha(fp), bound_source_hashes=d['bound_source_hashes'],
                         exact_process_journals_checked=d['exact_process_journals_checked'])
        gates[path] = entry
    root = Path(plan['output'])
    result = dict(status='verified_original_full_retained_timing_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), source_plan=str(pp), source_plan_sha256=sha(pp),
        frozen_pins_checked=len(plan['pins']), original_handles=handles, completion_gates=gates,
        output_root_present=root.exists(), completed_cohort_receipts=len(list((root / 'cohorts').glob('*.receipt.json'))),
        expected_cohorts=4340, full_candidate_census=5208000,
        producer_receipt_present=(root / 'receipt.json').exists(),
        independent_readback_present=(root / 'readback.json').exists(),
        production_fitting_launched=False, all_eight_aims_incomplete=True, gpu=False,
        scope='Original PID/create/command or invocation-specific terminal evidence, frozen plan/pins and '
              'current cgroup caps checked. Compact closure/archive presence and cohort counts only; '
              'no full source/artifact or timing replay, accepted fit, calibrated finish ETA or biological '
              'acceptance inferred. No job restart, limit change or original output mutation.')
    with args.output.open('x') as f:
        json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'original_handles'}, indent=2))


if __name__ == '__main__':
    main()
