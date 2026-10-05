#!/usr/bin/env python3
"""Run the unchanged full-proteome cataloger and bind its completed outputs."""
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
    args = parser.parse_args()
    assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text())
    pins = dict(plan['pins'])
    verify(pins)
    assert plan['expected_taxa'] == 526 and plan['expected_proteins'] == 5815847
    subprocess.run([sys.executable, 'scripts/catalog_whole_proteome_structures.py',
                    '--plan', str(args.plan)], check=True)
    root = Path(plan['output'])
    raw = root/'receipt.json'
    result = json.loads(raw.read_text())
    assert result['status'] == 'complete_whole_representative_proteome_exact_sequence_catalog'
    assert result['plan_sha256'] == sha(args.plan)
    assert result['taxa'] == plan['expected_taxa']
    assert result['proteins_screened'] == plan['expected_proteins']
    assert result['proteins_linked'] > 0 and result['unique_models'] > 0
    for name, digest in result['artifacts'].items():
        bind(pins, root/name, digest)
    for path in [args.plan, raw, Path(__file__)]:
        bind(pins, path)
    verify(pins)
    proof = dict(status='completed_full_afdb_catalog_pending_independent_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), raw_catalog_receipt=str(raw),
        raw_catalog_receipt_sha256=sha(raw), source_plan_sha256=sha(args.plan),
        taxa=result['taxa'], proteins_screened=result['proteins_screened'],
        proteins_linked=result['proteins_linked'], unique_models=result['unique_models'],
        taxa_with_models=result['taxa_with_models'],
        coordinate_bytes_verified=result['coordinate_bytes_verified'],
        inventory_records=result['inventory_records'], source_hashes=pins,
        scientific_eligibility=False, new_downloads=0, new_predictions=0, gpu=False,
        scope='Unchanged qualified cataloger checks every representative protein, latest accession '
              'status, exact sequence/length linkage, source-specific model ranking and selected '
              'coordinate byte hashes. All 526 taxa and 5,815,847 proteins retained. The raw '
              'producer receipt is immutable. Independent full link/model reconstruction and '
              'original execution closure remain required; no new coordinate-content, confidence '
              'or PAE qualification, structural clustering or biological acceptance.')
    with args.receipt.open('x') as handle:
        json.dump(proof, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in proof.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
