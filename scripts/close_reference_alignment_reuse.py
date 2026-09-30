#!/usr/bin/env python3
"""Archive full reuse proof/journal closure outside Git and publish its small locator."""
import argparse
import json
import subprocess
import sys
from pathlib import Path

from screen_duplication_alignment_reuse import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    def verify():
        for path, value in bindings.items():
            if sha(path) != value:
                raise ValueError('Changed control/source pin: ' + path)
    verify()
    source_plan = json.loads(Path(plan['source_plan']).read_text())
    root = Path(source_plan['output'])
    rp, ap = root / 'receipt.json', Path(plan['readback'])
    source, proof = json.loads(rp.read_text()), json.loads(ap.read_text())
    if (source['status'] != 'complete_full_reference_reuse_qualification_pending_independent_readback'
            or source['plan_sha256'] != sha(plan['source_plan'])
            or proof['status'] != 'passed_full_reference_input_checkpoint_numeric_reuse_readback'
            or proof['producer_receipt_sha256'] != sha(rp)):
        raise ValueError('Incomplete actual-data reuse qualification')
    fields = ['full_reference_pairs', 'directed_dispositions', 'counts', 'numerical_counts',
              'selected_source_dispositions', 'unique_model_mask_source_checks']
    summary = {k: source[k] for k in fields}
    if any(proof[k] != value for k, value in summary.items()):
        raise ValueError('Independent proof scope differs')
    if (source['full_reference_pairs'] != plan['full_reference_pairs']
            or source['directed_dispositions'] != 4 * source['full_reference_pairs']
            or sum(source['counts'].values()) != source['directed_dispositions']
            or source['counts'].get('new_native_measurement_pending') != plan['new_native_dispositions']):
        raise ValueError('Incomplete full design or altered source partition')
    archive_plan = root / 'completion_closure_plan.json'
    archive = root / 'completion_closure.json'
    inner = dict(output=str(archive), completed_status='complete_verified_full_reference_reuse_proof_archive',
                 evidence={'producer': dict(path=str(rp), expected=dict(status=source['status'], **summary)),
                           'reader': dict(path=str(ap), expected=dict(status=proof['status'], **summary))},
                 links=[{'from': 'reader', 'field': 'producer_receipt_sha256', 'to': str(rp)}],
                 launches=plan['launches'], pins={**bindings, str(rp): sha(rp), str(ap): sha(ap)},
                 summary=summary, scope=plan['scope'])
    with archive_plan.open('x') as handle:
        handle.write(json.dumps(inner, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py',
                    '--plan', str(archive_plan)], check=True)
    closure = json.loads(archive.read_text())
    if closure['status'] != inner['completed_status'] or closure['summary'] != summary or len(closure['services']) != 2:
        raise ValueError('Archived closure differs')
    verify()
    result = dict(status='complete_verified_full_reference_input_checkpoint_numeric_reuse',
                  **summary, full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive),
                  bound_source_and_artifact_hashes=len(closure['source_hashes']), exact_process_journals_checked=2,
                  producer_receipt=str(rp), producer_receipt_sha256=sha(rp),
                  independent_readback=str(ap), independent_readback_sha256=sha(ap),
                  source_plan=plan['source_plan'], source_plan_sha256=sha(plan['source_plan']),
                  completion_plan_sha256=sha(args.plan), scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['output']).open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
