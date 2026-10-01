#!/usr/bin/env python3
"""Close full new-background native, numeric/geometry and independent reader stages."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from expanded_background_measurement_sources import load_native
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args(); plan = json.loads(args.plan.read_text())
    _, root, native, _, _, _, bindings, _ = load_native(plan, include_positions=False); bind(bindings, args.plan)
    ap = Path(plan['assessment_plan']); qp = Path(plan['readback_plan'])
    assessment_plan, readback_plan = [json.loads(p.read_text()) for p in [ap, qp]]
    assert assessment_plan['native_plan'] == readback_plan['native_plan'] == plan['native_plan'] and readback_plan['assessment_plan'] == str(ap)
    dp, rp = Path(assessment_plan['output']) / 'receipt.json', Path(readback_plan['output'])
    diagnostic, reader = [json.loads(p.read_text()) for p in [dp, rp]]
    assert diagnostic['status'] == 'complete_expanded_background_numeric_geometry_pending_independent_readback' and diagnostic['plan_sha256'] == sha(ap)
    assert reader['status'] == 'passed_full_expanded_background_measurement_quaternion_readback' and reader['plan_sha256'] == sha(qp)
    assert diagnostic['native_receipt_sha256'] == sha(root / 'receipt.json') and reader['producer_receipt_sha256'] == sha(dp)
    fields = ['directed_dispositions', 'numerically_checked_alignments', 'counts', 'rmsd_status_counts', 'geometry_counts', 'numerical_eligibility_counts', 'maximum_rmsd_rounding_error', 'full_background_pairs', 'existing_catalog_pairs_pending_reuse']
    summary = {key: diagnostic[key] for key in fields}; assert all(reader[key] == value for key, value in summary.items())
    assert summary['directed_dispositions'] == native['directed_dispositions'] == 4 * plan['expected']['new_pairs']
    assert summary['counts'] == native['counts'] and summary['full_background_pairs'] == plan['expected']['full_pairs']
    assert summary['existing_catalog_pairs_pending_reuse'] == plan['expected']['pending_catalog_reuse_pairs']
    for path in [ap, qp, dp, rp]: bind(bindings, path)
    for proof in [diagnostic, reader]:
        assert proof['scientific_eligibility'] is False
        for path, digest in proof['source_hashes'].items(): bind(bindings, path, digest)
    verify(bindings)
    archive, inner_plan = root / 'measurement_completion_archive.json', root / 'measurement_completion_closure_plan.json'
    evidence = dict(native=dict(path=str(root / 'receipt.json'), expected=dict(status=native['status'], directed_dispositions=summary['directed_dispositions'])),
                    assessment=dict(path=str(dp), expected=dict(status=diagnostic['status'], **summary)), reader=dict(path=str(rp), expected=dict(status=reader['status'], **summary)))
    inner = dict(output=str(archive), completed_status='complete_verified_full_new_background_measurement_archive', evidence=evidence,
                 links=[{'from': 'assessment', 'field': 'native_receipt_sha256', 'to': str(root / 'receipt.json')}, {'from': 'reader', 'field': 'producer_receipt_sha256', 'to': str(dp)}],
                 pins=bindings, launches=plan['launches'], summary=summary, scope=plan['scope'])
    with inner_plan.open('x') as handle: handle.write(json.dumps(inner, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(inner_plan)], check=True)
    closed = json.loads(archive.read_text()); assert len(closed['services']) == 3
    result = dict(status='complete_verified_full_new_background_native_numeric_geometry', **summary,
                  new_pairs=plan['expected']['new_pairs'], full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive),
                  native_receipt=str(root / 'receipt.json'), native_receipt_sha256=sha(root / 'receipt.json'),
                  assessment_receipt=str(dp), assessment_receipt_sha256=sha(dp), independent_readback=str(rp), independent_readback_sha256=sha(rp),
                  bound_source_hashes=len(closed['source_hashes']), exact_process_journals_checked=3, completion_plan_sha256=sha(args.plan),
                  scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['output']).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__': main()
