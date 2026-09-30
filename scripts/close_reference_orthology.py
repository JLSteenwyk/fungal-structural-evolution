#!/usr/bin/env python3
"""Close the complete native reference membership export with both exact journals."""
import argparse
import json
import subprocess
import sys
from pathlib import Path

from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    source_plan = json.loads(Path(plan['source_plan']).read_text())
    root = Path(source_plan['output'])
    receipt_path, proof_path = root / 'receipt.json', Path(plan['readback'])
    source, proof = json.loads(receipt_path.read_text()), json.loads(proof_path.read_text())
    assert source['status'] == 'complete_full_reference_native_orthology_pending_readback'
    assert proof['status'] == 'passed_full_reference_native_orthology_readback'
    assert source['plan_sha256'] == proof['plan_sha256'] == sha(plan['source_plan'])
    assert proof['producer_receipt_sha256'] == sha(receipt_path)
    fields = ['target_contexts', 'unique_gene_pair_queries', 'reference_tie_records',
              'duplicate_reference_links', 'guides', 'context_reference_native_summary']
    summary = {k: source[k] for k in fields}
    assert all(proof[k] == v for k, v in summary.items())
    assert source['target_contexts'] == plan['target_contexts']
    assert source['reference_tie_records'] == plan['reference_tie_records']
    assert source['duplicate_reference_links'] == 2 * source['reference_tie_records']
    assert source['guide_contexts'] == {g: r['contexts'] for g, r in source_plan['resources']['guide_preflight'].items()}
    summary['guide_contexts'] = source['guide_contexts']
    inner = dict(output=plan['output'], completed_status='complete_verified_full_reference_native_orthology',
                 evidence=dict(producer=dict(path=str(receipt_path), expected=dict(status=source['status'], **summary)),
                               reader=dict(path=str(proof_path), expected=dict(status=proof['status'], **{k: v for k, v in summary.items() if k != 'guide_contexts'}))),
                 links=[{'from': 'reader', 'field': 'producer_receipt_sha256', 'to': str(receipt_path)}],
                 launches=plan['launches'], pins={**bindings, str(receipt_path): sha(receipt_path), str(proof_path): sha(proof_path)},
                 summary=summary, scope=plan['scope'])
    inner_plan = root / 'completion_closure_plan.json'
    with inner_plan.open('x') as handle:
        handle.write(json.dumps(inner, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(inner_plan)], check=True)
    for path, digest in bindings.items():
        assert sha(path) == digest, path


if __name__ == '__main__':
    main()
