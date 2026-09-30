#!/usr/bin/env python3
"""Close full reference order/coverage checks with both original process journals."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    for path, digest in bindings.items(): assert sha(path) == digest, path
    source_plan = json.loads(Path(plan['source_plan']).read_text()); root = Path(source_plan['output'])
    rp = root / 'receipt.json'; ap = Path(plan['readback']); r = json.loads(rp.read_text()); a = json.loads(ap.read_text())
    assert r['status'] == 'complete_full_reference_order_and_original_coverage_pending_independent_readback'
    assert a['status'] == 'passed_full_reference_order_and_original_coverage_readback'
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['source_plan']) and a['producer_receipt_sha256'] == sha(rp)
    fields = ['full_pairs', 'directed_dispositions', 'pair_mask_rows', 'pair_screen_decisions', 'order_summary_counts',
              'maximum_order_differences', 'pair_pass_counts', 'pair_exclusion_counts', 'pair_both_masks_pass_counts',
              'source_dispositions', 'native_status_counts', 'numerical_counts']
    summary = {k: r[k] for k in fields}; assert all(a[k] == v for k, v in summary.items())
    assert r['full_pairs'] == plan['full_pairs'] and r['directed_dispositions'] == 4 * plan['full_pairs']
    assert r['pair_mask_rows'] == 2 * plan['full_pairs'] and r['pair_screen_decisions'] == 12 * plan['full_pairs']
    assert sum(r['order_summary_counts'].values()) == r['pair_mask_rows']
    archive = root / 'completion_closure.json'; inner_plan = root / 'completion_closure_plan.json'
    inner = dict(output=str(archive), completed_status='complete_verified_full_reference_order_and_original_coverage_archive',
                 evidence=dict(producer=dict(path=str(rp), expected=dict(status=r['status'], **summary)), reader=dict(path=str(ap), expected=dict(status=a['status'], **summary))),
                 links=[{'from': 'reader', 'field': 'producer_receipt_sha256', 'to': str(rp)}], launches=plan['launches'],
                 pins={**bindings, str(rp): sha(rp), str(ap): sha(ap)}, summary=summary, scope=plan['scope'])
    with inner_plan.open('x') as handle: handle.write(json.dumps(inner, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(inner_plan)], check=True)
    closed = json.loads(archive.read_text()); assert closed['summary'] == summary and len(closed['services']) == 2
    result = dict(status='complete_verified_full_reference_order_and_original_coverage', **summary,
                  full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive), bound_source_and_artifact_hashes=len(closed['source_hashes']),
                  exact_process_journals_checked=2, producer_receipt=str(rp), producer_receipt_sha256=sha(rp),
                  independent_readback=str(ap), independent_readback_sha256=sha(ap), completion_plan_sha256=sha(args.plan),
                  scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['output']).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__': main()
