#!/usr/bin/env python3
"""Close full native-context coverage/assignment diagnostics with both process journals."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from reference_context_coverage_sources import POLICIES
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text()); bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    for path, h in bindings.items(): assert sha(path) == h, path
    config = json.loads(Path(plan['source_plan']).read_text()); root = Path(config['output']); rp = root / 'receipt.json'; ap = Path(plan['readback'])
    r, a = [json.loads(path.read_text()) for path in [rp, ap]]
    assert r['status'] == 'complete_full_reference_context_coverage_pending_independent_readback'
    assert a['status'] == 'passed_full_reference_context_coverage_sql_readback'
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['source_plan']) and a['producer_receipt_sha256'] == sha(rp)
    assert r['context_policy_flag_order'] == a['context_policy_flag_order'] == POLICIES
    fields = ['target_contexts', 'context_design_records', 'context_design_mask_records', 'reference_tie_records', 'duplicate_reference_links',
              'availability_side_links', 'context_screen_rows', 'context_policy_decisions', 'side_screen_decisions', 'summary_rows', 'guide_contexts', 'side_mask_work_counts']
    summary = {k: r[k] for k in fields}; assert all(a[k] == v for k, v in summary.items())
    for field in ['target_contexts', 'reference_tie_records', 'duplicate_reference_links', 'availability_side_links']: assert r[field] == config['expected'][field]
    n = r['target_contexts']; assert r['context_design_records'] == 2 * n and r['context_design_mask_records'] == 4 * n
    assert r['context_screen_rows'] == 24 * n and r['context_policy_decisions'] == 24 * n * len(POLICIES) and r['side_screen_decisions'] == 12 * r['duplicate_reference_links']
    assert r['summary_rows'] == 4 * 2 * 6 * len(POLICIES)
    inner = dict(output=plan['output'], completed_status='complete_verified_full_reference_context_coverage',
                 evidence=dict(producer=dict(path=str(rp), expected=dict(status=r['status'], **summary)), reader=dict(path=str(ap), expected=dict(status=a['status'], **summary))),
                 links=[{'from': 'reader', 'field': 'producer_receipt_sha256', 'to': str(rp)}], launches=plan['launches'],
                 pins={**bindings, str(rp): sha(rp), str(ap): sha(ap)}, summary=summary, scope=plan['scope'])
    inner_plan = root / 'completion_closure_plan.json'
    with inner_plan.open('x') as handle: handle.write(json.dumps(inner, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(inner_plan)], check=True)


if __name__ == '__main__': main()
