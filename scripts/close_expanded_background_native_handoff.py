#!/usr/bin/env python3
"""Close the full input/native-work handoff with complete hashes and original journals."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from expanded_background_input_sources import SUMMARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args(); plan = json.loads(args.plan.read_text()); bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    for path, digest in bindings.items(): assert sha(path) == digest, path
    config = json.loads(Path(plan['source_plan']).read_text()); root = Path(config['output']); rp, ap = root / 'receipt.json', root / 'readback.json'
    r, a = [json.loads(p.read_text()) for p in [rp, ap]]
    assert r['status'] == 'complete_full_expanded_background_native_input_handoff_pending_independent_readback' and a['status'] == 'passed_full_expanded_background_native_input_handoff_sql_readback'
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['source_plan']) and a['producer_receipt_sha256'] == sha(rp)
    summary = {key: r[key] for key in SUMMARY_FIELDS}; assert all(a[k] == v for k, v in summary.items())
    for key in ['active_models', 'full_pairs', 'new_pairs', 'pending_catalog_reuse_pairs']: assert summary[key] == config['expected'][key]
    archive = root / 'completion_archive.json'; inner_plan = root / 'completion_closure_plan.json'
    inner = dict(output=str(archive), completed_status='complete_verified_full_expanded_background_native_input_handoff_archive',
                 evidence=dict(producer=dict(path=str(rp), expected=dict(status=r['status'], **summary)), reader=dict(path=str(ap), expected=dict(status=a['status'], **summary))),
                 links=[{'from': 'reader', 'field': 'producer_receipt_sha256', 'to': str(rp)}], launches=plan['launches'], pins={**bindings, str(rp): sha(rp), str(ap): sha(ap)}, summary=summary, scope=plan['scope'])
    with inner_plan.open('x') as handle: handle.write(json.dumps(inner, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(inner_plan)], check=True)
    closed = json.loads(archive.read_text()); assert len(closed['services']) == 2
    result = dict(status='complete_verified_full_expanded_background_native_input_handoff', **summary, producer_receipt=str(rp), producer_receipt_sha256=sha(rp),
                  independent_readback=str(ap), independent_readback_sha256=sha(ap), full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive),
                  bound_source_hashes=len(closed['source_hashes']), exact_process_journals_checked=2, scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['output']).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__': main()
