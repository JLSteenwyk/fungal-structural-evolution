#!/usr/bin/env python3
"""Verify and close all four new-reference measurement stages before full union."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from reference_measurement_union_sources import native_bundle, verify, bind
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    native = native_bundle(plan, bindings)
    verify(bindings)
    nr, dr, gr = [native[k] for k in ['native_receipt', 'diagnostic_receipt', 'geometry_receipt']]
    summary = dict(full_pairs=plan['full_pairs'], new_pairs=plan['new_pairs'],
                   directed_dispositions=nr['directed_dispositions'], native_counts=nr['counts'],
                   numerically_checked_alignments=dr['numerically_checked_alignments'],
                   rmsd_status_counts=dr['rmsd_status_counts'], geometry_counts=gr['counts'])
    root = native['root']
    dp = Path(json.loads(Path(plan['native_diagnostic_plan']).read_text())['output'])
    gp = Path(json.loads(Path(plan['native_geometry_plan']).read_text())['output'])
    evidence = dict(native=dict(path=str(root/'receipt.json'), expected=dict(status=nr['status'], directed_dispositions=summary['directed_dispositions'])),
                    diagnostic=dict(path=str(dp/'receipt.json'), expected=dict(status=dr['status'], directed_dispositions=summary['directed_dispositions'])),
                    geometry=dict(path=str(gp/'receipt.json'), expected=dict(status=gr['status'], alignments=summary['numerically_checked_alignments'])),
                    reader=dict(path=plan['native_readback'], expected=dict(status='passed_full_reference_geometry_readback', alignments_checked=summary['numerically_checked_alignments'])))
    archive = root / 'measurement_completion_closure.json'
    inner = dict(output=str(archive), completed_status='complete_verified_new_reference_native_measurement_archive',
                 evidence=evidence, pins=bindings, launches=plan['launches'], summary=summary, scope=plan['scope'],
                 links=[{'from':'diagnostic','field':'producer_receipt_sha256','to':str(root/'receipt.json')},
                        {'from':'geometry','field':'diagnostic_receipt_sha256','to':str(dp/'receipt.json')},
                        {'from':'reader','field':'producer_receipt_sha256','to':str(gp/'receipt.json')}])
    for spec in evidence.values(): bind(inner['pins'], spec['path'])
    inner_plan = root / 'measurement_completion_closure_plan.json'
    with inner_plan.open('x') as handle: handle.write(json.dumps(inner, indent=2)+'\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(inner_plan)], check=True)
    closed = json.loads(archive.read_text())
    assert len(closed['services']) == 4 and closed['summary'] == summary
    result = dict(status='complete_verified_new_reference_native_alignment_numeric_geometry', **summary,
                  full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive),
                  bound_source_and_artifact_hashes=len(closed['source_hashes']), exact_process_journals_checked=4,
                  completion_plan_sha256=sha(args.plan), scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['output']).open('x') as handle: handle.write(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__': main()
