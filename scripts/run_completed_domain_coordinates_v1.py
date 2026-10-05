#!/usr/bin/env python3
"""Extract the entire refreshed domain manifest after its original readback closes."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

import Bio

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text())
    assert '.'.join(map(str, sys.version_info[:3])) == plan['environment']['python']
    assert Bio.__version__ == plan['environment']['biopython']
    pins = dict(plan['pins'])
    verify(pins)
    audit_path = Path(plan['manifest_readback_receipt'])
    audit = json.loads(audit_path.read_text())
    assert audit['status'] == 'passed_full_domain_interval_manifest_sql_union_readback'
    assert audit['models'] == plan['expected_models']
    assert audit['intervals'] == plan['expected_intervals']
    assert audit['original_producer_tool_terminal_exit_code'] == 1
    assert audit['original_producer_native_terminal_exit_code'] is None
    for kind in ('reporting_review', 'readback'):
        tp = Path(plan['manifest_'+kind+'_transport'])
        t = json.loads(tp.read_text())
        vp = Path(plan['manifest_'+kind+'_receipt'])
        assert t['status'] == 'verified_original_software_wait_and_whole_wrapper_payloads'
        assert t['validation_sha256'] == sha(vp)
        assert t['original_tool_terminal_exit_code'] == 0
        assert t['whole_wrapper_initial_and_terminal_payloads_matched']
        assert t['manager_start_records'] == t['manager_completion_records'] == 1
        for path, digest in t['source_hashes'].items():
            bind(pins, path, digest)
    verify(pins)
    subprocess.run([sys.executable, 'scripts/extract_domain_coordinates.py',
                    '--plan', str(args.plan)], check=True)
    root = Path(plan['output'])
    raw_path = root/'receipt.json'
    raw = json.loads(raw_path.read_text())
    assert raw['status'] == 'complete_domain_coordinate_dispositions_pending_independent_readback'
    assert raw['plan_sha256'] == sha(args.plan)
    assert raw['counts']['models'] == plan['expected_models']
    assert raw['counts']['intervals'] == plan['expected_intervals']
    assert raw['counts']['exported']+raw['counts']['rejected'] == plan['expected_intervals']
    for proof in raw['proofs']:
        bind(pins, proof['job'], proof['receipt']['job_sha256'])
        for name, digest in proof['receipt']['artifacts'].items():
            bind(pins, root/'shards'/name, digest)
        bind(pins, root/'shards'/(Path(proof['job']).stem+'.receipt.json'))
    for path in (args.plan, raw_path, audit_path, Path(__file__)):
        bind(pins, path)
    verify(pins)
    result = dict(status='completed_full_refreshed_domain_coordinates_pending_independent_atom_readback',
                  checked_utc=datetime.now(timezone.utc).isoformat(), counts=raw['counts'],
                  shards=raw['shards'], raw_receipt_sha256=sha(raw_path),
                  source_hashes=pins, scientific_eligibility=False, gpu=False, new_predictions=0,
                  scope='All refreshed eligible domain intervals attempted with unchanged source, '
                        'atom, sequence, serialization and confidence checks. Every rejection and '
                        'missing-backbone disposition retained. Original reporting failure preserved; '
                        'both actual original review/readback waits closed. No native rerun or assumed exit. '
                        'Full independent archive/atom readback remains required; no PAE, biological '
                        'boundary, homology, function or evolutionary-event acceptance.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
