#!/usr/bin/env python3
"""Close full original-residue maps/readback in an outside-Git hash archive."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from readback_full_triad_common_residues import SUMMARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text()); bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    for path, h in bindings.items(): assert sha(path) == h, path
    config = json.loads(Path(plan['source_plan']).read_text()); root = Path(config['output']); rp = root / 'receipt.json'; ap = root / 'readback.json'
    r, a = [json.loads(path.read_text()) for path in [rp, ap]]
    assert r['status'] == 'complete_full_triad_original_residue_mapping_pending_independent_readback'
    assert a['status'] == 'passed_full_triad_original_residue_mapping_readback'
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['source_plan']) and a['producer_receipt_sha256'] == sha(rp)
    summary = {k: r[k] for k in SUMMARY_FIELDS}; assert all(a[k] == v for k, v in summary.items())
    assert summary['mask_order_states'] == config['expected']['potential_correspondence_states'] == 16 * config['expected']['correspondence_work_triads']
    for field in ['target_contexts', 'reference_tie_records', 'duplicate_reference_links', 'unique_ordered_model_triads', 'correspondence_work_triads', 'all_input_dispositions']:
        assert summary[field] == config['expected'][field]
    assert summary['cycle_consistent_residue_occurrences'] <= summary['common_reference_residue_occurrences']
    archive_path = root / 'completion_archive.json'
    inner = dict(output=str(archive_path), completed_status='complete_verified_full_triad_original_residue_mapping_archive',
                 evidence=dict(producer=dict(path=str(rp), expected=dict(status=r['status'], **summary)), reader=dict(path=str(ap), expected=dict(status=a['status'], **summary))),
                 links=[{'from': 'reader', 'field': 'producer_receipt_sha256', 'to': str(rp)}], launches=plan['launches'],
                 pins={**bindings, str(rp): sha(rp), str(ap): sha(ap)}, summary=summary, scope=plan['scope'])
    inner_plan = root / 'completion_closure_plan.json'
    with inner_plan.open('x') as handle: handle.write(json.dumps(inner, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(inner_plan)], check=True)
    archive = json.loads(archive_path.read_text()); assert len(archive['services']) == 2
    result = dict(status='complete_verified_full_triad_original_residue_mapping', **summary,
                  producer_receipt=str(rp), producer_receipt_sha256=sha(rp), independent_readback=str(ap), independent_readback_sha256=sha(ap),
                  full_hash_archive=str(archive_path), full_hash_archive_sha256=sha(archive_path), bound_source_hashes=len(archive['source_hashes']),
                  exact_process_journals_checked=2, scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['output']).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__': main()
