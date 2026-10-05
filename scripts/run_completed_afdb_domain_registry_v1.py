#!/usr/bin/env python3
"""Run unchanged full domain registry only from a closed complete catalog."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text()); pins = dict(plan['pins']); verify(pins)
    closure_path = Path(plan['catalog_closure'])
    closure = json.loads(closure_path.read_text())
    assert closure['status'] == 'complete_verified_full_completed_retrieval_catalog_refresh'
    assert closure['taxa'] == 526 and closure['proteins_screened'] == 5815847
    assert closure['unique_models'] == plan['expected_models']
    assert closure['proteins_linked'] == plan['expected_protein_links']
    comparison = json.loads(Path(plan['catalog_comparison_closure']).read_text())
    assert comparison['status'] == 'complete_verified_full_afdb_catalog_comparison'
    assert comparison['new_links'] == plan['expected_protein_links']
    for path in [plan['models'], plan['protein_links']]:
        assert pins[path] == closure['source_hashes'][path]
    subprocess.run([sys.executable, 'scripts/build_structure_domain_registry.py', '--plan', str(args.plan)], check=True)
    root = Path(plan['output']); raw_path = root / 'receipt.json'
    raw = json.loads(raw_path.read_text())
    assert raw['status'] == 'complete_structure_domain_interval_registry_pending_independent_readback'
    assert raw['plan_sha256'] == sha(args.plan)
    assert raw['counts']['models'] == plan['expected_models']
    assert raw['counts']['protein_links'] == plan['expected_protein_links']
    bind(pins, root / 'structure_domains.sqlite', raw['database_sha256'])
    for path in [args.plan, raw_path, Path(__file__), root / 'state.json']: bind(pins, path)
    verify(pins)
    result = dict(status='completed_full_refreshed_afdb_domain_registry_pending_independent_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), counts=raw['counts'],
        candidate_domain_intervals_by_policy=raw['candidate_domain_intervals_by_policy'],
        original_taxa=526, original_representative_proteins=5815847,
        source_unlinked_proteins=closure['proteins_without_catalog_model'],
        raw_registry_receipt_sha256=sha(raw_path), database_sha256=raw['database_sha256'],
        source_hashes=pins, scientific_eligibility=False, gpu=False, new_predictions=0,
        scope='Unchanged qualified whole catalog model/annotation interval join over every selected model '
              'and all protein links. Existing four architecture policy alternatives and HMM/type/length '
              'rules retained. Catalog missingness remains explicit outside model registry. Full independent '
              'interval/protein reconstruction and original actual execution closure required. No residue '
              'confidence/PAE, structural domain validation, homology or evolutionary event acceptance.')
    with args.receipt.open('x') as handle: json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__': main()
