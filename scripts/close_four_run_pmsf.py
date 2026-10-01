#!/usr/bin/env python3
"""Close full crossed PMSF sensitivity with complete grids and original journals."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from four_run_pmsf_sources import SUMMARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text()); config = json.loads(Path(plan['source_plan']).read_text())
    root = Path(config['output']); rp = root / 'receipt.json'; ap = root / 'readback.json'; r, a = [json.loads(p.read_text()) for p in [rp, ap]]
    assert r['status'] == 'complete_full_four_run_pmsf_ML_and_consensus_sensitivity_pending_readback'
    assert a['status'] == 'passed_full_four_run_pmsf_ML_consensus_split_support_readback'
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['source_plan']) and a['producer_receipt_sha256'] == sha(rp)
    summary = {key: r[key] for key in SUMMARY_FIELDS}; assert all(a[key] == value for key, value in summary.items())
    assert all(summary[key] == config['expected'][key] for key in ['taxa', 'internal_splits_per_view', 'tree_views', 'comparison_rows'])
    bindings = {str(args.plan): sha(args.plan), **plan['pins'], str(rp): sha(rp), str(ap): sha(ap)}
    inner = dict(output=plan['output'], completed_status='complete_verified_full_four_run_pmsf_ML_consensus_sensitivity',
                 evidence=dict(producer=dict(path=str(rp), expected=dict(status=r['status'], **summary)), reader=dict(path=str(ap), expected=dict(status=a['status'], **summary))),
                 links=[{'from': 'reader', 'field': 'producer_receipt_sha256', 'to': str(rp)}], launches=plan['launches'],
                 pins=bindings, summary=summary, scope=config['scope'])
    cp = root / 'completion_closure_plan.json'
    with cp.open('x') as f: f.write(json.dumps(inner, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(cp)], check=True)


if __name__ == '__main__': main()
