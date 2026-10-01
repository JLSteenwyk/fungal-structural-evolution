#!/usr/bin/env python3
"""Close full identical-residue geometry with hashes and exact original journals."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from full_triad_fit_sources import SUMMARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text()); bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    for path, digest in bindings.items(): assert sha(path) == digest, path
    config = json.loads(Path(plan['source_plan']).read_text()); root = Path(config['output']); rp = root / 'receipt.json'; ap = root / 'readback.json'
    r, a = [json.loads(path.read_text()) for path in [rp, ap]]
    assert r['status'] == 'complete_full_triad_same_residue_fits_pending_independent_readback'
    assert a['status'] == 'passed_full_triad_same_residue_quaternion_readback'
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['source_plan']) and a['producer_receipt_sha256'] == sha(rp)
    summary = {k: r[k] for k in SUMMARY_FIELDS}; assert all(a[k] == value for k, value in summary.items())
    assert summary['mapping_states'] == config['expected']['mask_order_states'] and summary['fit_rows'] == 2 * summary['mapping_states']
    archive_path = root / 'completion_archive.json'
    inner = dict(output=str(archive_path), completed_status='complete_verified_full_triad_same_residue_fit_archive',
                 evidence=dict(producer=dict(path=str(rp), expected=dict(status=r['status'], **summary)), reader=dict(path=str(ap), expected=dict(status=a['status'], **summary))),
                 links=[{'from': 'reader', 'field': 'producer_receipt_sha256', 'to': str(rp)}], launches=plan['launches'],
                 pins={**bindings, str(rp): sha(rp), str(ap): sha(ap)}, summary=summary, scope=plan['scope'])
    inner_plan = root / 'completion_closure_plan.json'
    with inner_plan.open('x') as f: f.write(json.dumps(inner, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(inner_plan)], check=True)
    archive = json.loads(archive_path.read_text()); assert len(archive['services']) == 2
    result = dict(status='complete_verified_full_triad_same_residue_geometry', **summary,
                  maximum_absolute_rmsd_or_contrast_difference=a['maximum_absolute_rmsd_or_contrast_difference'],
                  producer_receipt=str(rp), producer_receipt_sha256=sha(rp), independent_readback=str(ap), independent_readback_sha256=sha(ap),
                  full_hash_archive=str(archive_path), full_hash_archive_sha256=sha(archive_path), bound_source_hashes=len(archive['source_hashes']),
                  exact_process_journals_checked=2, scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['output']).open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__': main()
