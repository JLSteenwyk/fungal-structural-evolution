#!/usr/bin/env python3
"""Retranslate every original CDS under predeclared inherited and snapshot codes."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
from pathlib import Path
import time

from full_unmodified_cds_translation_v1 import audit_taxon
from reference_measurement_union_sources import bind, verify


FIELDS = ['target_role_status_counts', 'inherited_status_cross_counts', 'changed_codon_role_counts',
          'product_role_status_counts', 'coding_model_translation_cross_counts']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['plan', 'receipt']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text())
    assert plan['expected_taxa'] == len(plan['entries']) == 526
    verify(plan['pins'])
    root = Path(plan['output'])
    root.mkdir(exist_ok=False)
    start, reports = time.monotonic(), []
    pins = dict(plan['pins'])
    with ProcessPoolExecutor(max_workers=plan['cpu']) as pool:
        futures = [pool.submit(audit_taxon, entry, root) for entry in plan['entries']]
        for future in as_completed(futures):
            report = future.result()
            reports.append(report)
            for path in report['artifact_paths']:
                bind(pins, path)
            state = dict(stage='unmodified_original_cds_fixed_code_translation', completed_taxa=len(reports),
                expected_taxa=526, source_products=sum(r['source_products'] for r in reports),
                target_records=sum(r['target_records'] for r in reports),
                changed_codon_rows=sum(r['changed_codon_rows'] for r in reports),
                elapsed_seconds=time.monotonic() - start)
            temporary = root / 'state.partial'
            temporary.write_text(json.dumps(state, indent=2) + '\n')
            temporary.replace(root / 'state.json')
            print(json.dumps(state), flush=True)
    assert len(reports) == len({r['taxon_id'] for r in reports}) == 526
    totals = {field: sum(r[field] for r in reports) for field in ['source_products', 'selected_representatives',
               'target_records', 'original_dna_bases', 'changed_codon_rows']}
    for field in ['source_products', 'selected_representatives', 'target_records']:
        assert totals[field] == plan[field]
    aggregates = {field: Counter() for field in FIELDS}
    for report in reports:
        for field in FIELDS:
            aggregates[field].update(report[field])
    bind(pins, args.plan)
    verify(pins)
    result = dict(status='complete_full_unmodified_cds_fixed_code_translation_pending_independent_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, **totals,
        taxa_reports=sorted(reports, key=lambda r: r['taxon_id']),
        **{f: dict(v) for f, v in aggregates.items()}, elapsed_seconds=time.monotonic() - start,
        source_hashes=pins, scientific_eligibility=False, biological_codon_eligibility=False,
        genetic_code_admission=False, dna_modified=False, protein_modified=False, best_code_selected=False,
        gpu=False, new_predictions=0, all_eight_aims_incomplete=True,
        scope='All original targets independently translated from their unmodified FASTAs under fixed '
              'inherited and pinned nuclear/mitochondrial choices. All original source products, alternatives, '
              'unlinked targets and no-target cases retained. Nontriplet sequences are not translated; at most '
              'one terminal stop removed for comparison. No frame, initiator, recoding, phase or boundary '
              'repair, no derived target substitution, no best-code selection or compartment/selection inference. '
              'Complete independent codon-table replay and original transport closure still required.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'taxa_reports']}, indent=2))


if __name__ == '__main__':
    main()
