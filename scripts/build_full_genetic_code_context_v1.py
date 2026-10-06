#!/usr/bin/env python3
"""Compare all inherited coding test codes to frozen taxonomy source context."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
from pathlib import Path
import time

from genetic_code_context_v1 import audit_taxon
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    args = p.parse_args()
    assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text())
    verify(plan['pins'])
    root = Path(plan['output'])
    root.mkdir(exist_ok=False)
    start, reports = time.monotonic(), []
    pins = dict(plan['pins'])
    with ProcessPoolExecutor(max_workers=plan['cpu']) as pool:
        futures = [pool.submit(audit_taxon, e, root) for e in plan['entries']]
        for future in as_completed(futures):
            report = future.result()
            reports.append(report)
            for path in report['artifact_paths']:
                bind(pins, path)
            state = dict(completed_taxa=len(reports), expected_taxa=526,
                         source_products=sum(r['source_products'] for r in reports),
                         elapsed_seconds=time.monotonic() - start)
            temporary = root / 'state.partial'
            temporary.write_text(json.dumps(state, indent=2) + '\n')
            temporary.replace(root / 'state.json')
            print(json.dumps(state), flush=True)
    assert len(reports) == len({r['taxon_id'] for r in reports}) == 526
    totals = {f: sum(r[f] for r in reports) for f in ['source_products', 'selected_representatives', 'target_records']}
    assert totals == dict(source_products=5927745, selected_representatives=5815847, target_records=5923039)
    aggregates = {f: Counter() for f in ['target_relationship_counts', 'product_classification_counts',
                                        'coding_model_code_cross_counts']}
    for report in reports:
        for field in aggregates:
            aggregates[field].update(report[field])
    bind(pins, args.plan)
    verify(pins)
    result = dict(status='complete_full_genetic_code_context_pending_independent_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, **totals,
        taxa_reports=sorted(reports, key=lambda r: r['taxon_id']),
        **{f: dict(v) for f, v in aggregates.items()}, elapsed_seconds=time.monotonic() - start,
        source_hashes=pins, scientific_eligibility=False, genetic_code_admission=False,
        biological_codon_eligibility=False, translations_recomputed=False, gpu=False, new_predictions=0,
        all_eight_aims_incomplete=True,
        scope='Full original CDS target/product test-code relationships to pinned selected taxonomy. '
              'Every alternative, source review and derived evidence retained. No code correction, '
              'translation, compartment inference, current-taxonomy claim, codon alignment, selection '
              'or any admitted evolutionary result. Independent complete replay remains required.')
    with args.receipt.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'taxa_reports']}, indent=2))


if __name__ == '__main__':
    main()
