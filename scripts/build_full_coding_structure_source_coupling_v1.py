#!/usr/bin/env python3
"""Link every original CDS/product to inherited translation and model evidence."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
from pathlib import Path
import time

from coding_structure_source_coupling_v1 import couple_taxon
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
    start = time.monotonic()
    reports = []
    pins = dict(plan['pins'])
    with ProcessPoolExecutor(max_workers=plan['cpu']) as pool:
        futures = [pool.submit(couple_taxon, entry, plan['availability_database'], root) for entry in plan['entries']]
        for future in as_completed(futures):
            report = future.result()
            reports.append(report)
            for path in report['artifact_paths']:
                bind(pins, path)
            state = dict(stage='full_coding_structure_source_coupling', completed_taxa=len(reports), expected_taxa=526,
                         source_products=sum(r['source_products'] for r in reports), elapsed_seconds=time.monotonic() - start)
            temporary = root / 'state.partial'
            temporary.write_text(json.dumps(state, indent=2) + '\n')
            temporary.replace(root / 'state.json')
            print(json.dumps(state), flush=True)
    assert len(reports) == len({r['taxon_id'] for r in reports}) == 526
    assert sum(r['source_products'] for r in reports) == 5927745
    assert sum(r['selected_representatives'] for r in reports) == 5815847
    assert sum(r['target_records'] for r in reports) == 5923039
    products, targets, cross = Counter(), Counter(), Counter()
    for report in reports:
        products.update(report['product_status_counts'])
        targets.update(report['target_status_counts'])
        cross.update(report['coding_structure_cross_counts'])
    availability = Counter()
    for key, count in cross.items():
        availability[key.rsplit(':', 1)[1]] += count
    assert dict(availability) == dict(afdb_only=2994146, esmfold_only=24801, both=722,
                                     neither=2796178, alternative_not_assigned=111898)
    bind(pins, args.plan)
    verify(pins)
    result = dict(status='complete_full_coding_structure_source_coupling_pending_independent_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, source_products=5927745,
        selected_representatives=5815847, target_records=5923039,
        taxa_reports=sorted(reports, key=lambda r: r['taxon_id']),
        product_status_counts=dict(products), target_status_counts=dict(targets),
        coding_structure_cross_counts=dict(cross), representative_availability_and_alternative_counts=dict(availability),
        elapsed_seconds=time.monotonic() - start, source_hashes=pins,
        scientific_eligibility=False, biological_codon_eligibility=False,
        inherited_translation_independently_recomputed=False, gpu=False, new_predictions=0,
        all_eight_aims_incomplete=True,
        scope='All526entries/all5927745original products/all5923039original CDS targets source-bound to existing '
              'translation classifications and complete representative model availability. Separate projected/extracted '
              'codon evidence never replaces original target evidence. Code assumptions, gene/isoform uncertainties, '
              'all alternatives and review cases remain explicit. Translation is inherited and not independently '
              'recomputed here; source joins require full independent replay. No confidence, copy, homology, taxonomy, '
              'contamination, alignment/divergence, selection or evolutionary admission.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'taxa_reports']}, indent=2))


if __name__ == '__main__':
    main()
