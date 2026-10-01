#!/usr/bin/env python3
"""Close full current background reuse qualification and ten original source/check journals."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); args = parser.parse_args()
    plan = json.loads(args.plan.read_text()); config = json.loads(Path(plan['source_plan']).read_text()); root = Path(config['output']); rp, ap = root / 'receipt.json', Path(plan['readback'])
    r, a = [json.loads(path.read_text()) for path in [rp, ap]]; bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    assert r['status'] == 'complete_full_expanded_background_reuse_qualification_pending_independent_readback' and r['plan_sha256'] == sha(plan['source_plan'])
    assert a['status'] == 'passed_full_expanded_background_input_checkpoint_numeric_quaternion_reuse_readback' and a['producer_receipt_sha256'] == sha(rp) and a['plan_sha256'] == sha(plan['source_plan'])
    fields = ['full_background_pairs', 'directed_dispositions', 'counts', 'numerical_counts', 'selected_source_dispositions', 'unique_model_mask_source_checks', 'numerically_reconstructed_reuse_alignments']
    summary = {key: r[key] for key in fields}; assert all(a[key] == value for key, value in summary.items())
    assert summary['full_background_pairs'] == config['expected']['full_pairs'] and summary['directed_dispositions'] == 4 * summary['full_background_pairs']
    assert summary['counts'].get('new_native_measurement_pending') == 4 * config['expected']['new_pairs']
    assert set(summary['counts']) <= {'new_native_measurement_pending', 'incompatible_inputs_require_new_measurement', 'verified_identical_input_checkpoint_and_retained_disposition'}
    assert sum(summary['counts'].values()) == summary['directed_dispositions']
    evidence = dict(producer=dict(path=str(rp), expected=dict(status=r['status'], **summary)), reader=dict(path=str(ap), expected=dict(status=a['status'], **summary)))
    links = [{'from': 'reader', 'field': 'producer_receipt_sha256', 'to': str(rp)}]
    for label, spec in config['sources'].items():
        paths = [Path(json.loads(Path(spec['alignment_plan']).read_text())['output']) / 'receipt.json', Path(spec['diagnostic']) / 'receipt.json', Path(spec['geometry']) / 'receipt.json', Path(spec['geometry_readback'])]
        states = [spec['alignment_status'], spec['diagnostic_status'], spec['geometry_status'], spec['geometry_readback_status']]
        for index, (path, status) in enumerate(zip(paths, states)): evidence[label + str(index)] = dict(path=str(path), expected=dict(status=status)); bind(bindings, path)
        links.extend([{'from': label + '1', 'field': 'producer_receipt_sha256', 'to': str(paths[0])}, {'from': label + '2', 'field': 'diagnostic_receipt_sha256', 'to': str(paths[1])}, {'from': label + '3', 'field': 'producer_receipt_sha256', 'to': str(paths[2])}])
    for path in [rp, ap]: bind(bindings, path)
    for proof in [r, a]:
        assert proof['scientific_eligibility'] is False
        for path, digest in proof['source_hashes'].items(): bind(bindings, path, digest)
    verify(bindings)
    archive, inner_plan = root / 'completion_archive.json', root / 'completion_closure_plan.json'
    inner = dict(output=str(archive), completed_status='complete_verified_full_background_reuse_disposition_archive', evidence=evidence, links=links,
                 launches=plan['launches'], pins=bindings, summary=summary, scope=plan['scope'])
    with inner_plan.open('x') as handle: handle.write(json.dumps(inner, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(inner_plan)], check=True)
    closed = json.loads(archive.read_text()); assert len(closed['services']) == 10
    result = dict(status='complete_verified_full_background_input_checkpoint_reuse_with_dispositions', **summary, full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive),
        bound_source_hashes=len(closed['source_hashes']), exact_process_journals_checked=10, producer_receipt=str(rp), producer_receipt_sha256=sha(rp),
        independent_readback=str(ap), independent_readback_sha256=sha(ap), source_plan=plan['source_plan'], source_plan_sha256=sha(plan['source_plan']),
        completion_plan_sha256=sha(args.plan), scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['output']).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__': main()
