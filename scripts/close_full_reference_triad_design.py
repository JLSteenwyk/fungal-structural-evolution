#!/usr/bin/env python3
"""Close full three-edge work design with independent readback and exact journals."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from reference_triad_design_sources import SUMMARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text()); bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    for path, h in bindings.items(): assert sha(path) == h, path
    config = json.loads(Path(plan['source_plan']).read_text()); root = Path(config['output']); rp = root / 'receipt.json'; ap = Path(plan['readback'])
    r, a = [json.loads(path.read_text()) for path in [rp, ap]]
    assert r['status'] == 'complete_full_reference_triad_work_design_pending_independent_readback'
    assert a['status'] == 'passed_full_reference_triad_work_design_sql_readback'
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['source_plan']) and a['producer_receipt_sha256'] == sha(rp)
    summary = {k: r[k] for k in SUMMARY_FIELDS}; assert all(a[k] == v for k, v in summary.items())
    for field in ['target_contexts', 'reference_tie_records', 'duplicate_reference_links']: assert summary[field] == config['expected'][field]
    assert summary['availability_side_links_preserved'] == config['expected']['availability_side_links']
    assert summary['context_design_records'] == 2 * summary['target_contexts'] and summary['duplicate_reference_links'] == 2 * summary['reference_tie_records']
    assert summary['potential_correspondence_states'] == 16 * summary['correspondence_work_triads']
    inner = dict(output=plan['output'], completed_status='complete_verified_full_reference_triad_work_design',
                 evidence=dict(producer=dict(path=str(rp), expected=dict(status=r['status'], **summary)), reader=dict(path=str(ap), expected=dict(status=a['status'], **summary))),
                 links=[dict(**{'from': 'reader'}, field='producer_receipt_sha256', to=str(rp))], launches=plan['launches'],
                 pins={**bindings, str(rp): sha(rp), str(ap): sha(ap)}, summary=summary, scope=plan['scope'])
    inner_plan = root / 'completion_closure_plan.json'
    with inner_plan.open('x') as handle: handle.write(json.dumps(inner, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(inner_plan)], check=True)


if __name__ == '__main__': main()
