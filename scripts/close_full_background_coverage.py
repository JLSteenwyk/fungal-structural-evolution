#!/usr/bin/env python3
"""Require full coverage proof, all bound sources and both exact original journals."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from full_background_coverage_sources import load, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); args = parser.parse_args()
    plan = json.loads(args.plan.read_text()); source_plan = json.loads(Path(plan['source_plan']).read_text()); _, bindings = load(source_plan, plan['source_plan']); bind(bindings, args.plan)
    for path, digest in plan['pins'].items(): bind(bindings, path, digest)
    root = Path(source_plan['output']); rp, ap = root / 'receipt.json', Path(plan['readback']); r, a = [json.loads(p.read_text()) for p in [rp, ap]]
    assert r['status'] == 'complete_full_background_coverage_pending_independent_readback' and r['plan_sha256'] == sha(plan['source_plan'])
    assert a['status'] == 'passed_full_background_coverage_sql_decimal_readback' and a['plan_sha256'] == sha(plan['source_plan']) and a['producer_receipt_sha256'] == sha(rp)
    summary = {key: r[key] for key in SUMMARY_FIELDS}; assert all(a[key] == value for key, value in summary.items())
    assert summary['pairs'] == source_plan['expected']['pairs'] and summary['pair_mask_rows'] == 2 * summary['pairs'] and summary['directed_source_states'] == 4 * summary['pairs'] and summary['screen_decisions'] == 2 * summary['pairs'] * len(source_plan['screens'])
    for p in [rp, ap]: bind(bindings, p)
    for record in [r, a]:
        assert record['scientific_eligibility'] is False
        for path, digest in record['source_hashes'].items(): bind(bindings, path, digest)
    verify(bindings); archive, inner_plan = root / 'completion_archive.json', root / 'completion_closure_plan.json'
    inner = dict(output=str(archive), completed_status='complete_verified_full_background_coverage_archive', evidence=dict(producer=dict(path=str(rp), expected=dict(status=r['status'], **summary)), reader=dict(path=str(ap), expected=dict(status=a['status'], **summary))),
                 links=[{'from': 'reader', 'field': 'producer_receipt_sha256', 'to': str(rp)}], launches=plan['launches'], pins=bindings, summary=summary, scope=plan['scope'])
    with inner_plan.open('x') as handle: handle.write(json.dumps(inner, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(inner_plan)], check=True)
    closed = json.loads(archive.read_text()); assert len(closed['services']) == 2
    result = dict(status='complete_verified_full_background_original_length_coverage', **summary, full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive), bound_source_hashes=len(closed['source_hashes']), exact_process_journals_checked=2,
                  producer_receipt=str(rp), producer_receipt_sha256=sha(rp), independent_readback=str(ap), independent_readback_sha256=sha(ap), source_plan=plan['source_plan'], source_plan_sha256=sha(plan['source_plan']), completion_plan_sha256=sha(args.plan), scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['output']).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
