#!/usr/bin/env python3
"""Observe the original full catalog producer and queued independent reader."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    launches = ['metadata/completed_afdb_catalog_launch_20261005_v1.json',
                'metadata/completed_afdb_catalog_readback_launch_20261005_v1.json']
    limits = {'cpu.max': '200000 100000', 'memory.max': str(96*2**30),
              'memory.swap.max': '0'}
    handles = [observe(path, limits) for path in launches]
    pp = Path('metadata/completed_afdb_catalog_plan_20261005_v1.json')
    plan = json.loads(pp.read_text())
    root = Path(plan['output'])
    stdout = Path('metadata/completed_afdb_catalog_execution_20261005_v1/stdout.log')
    stderr = stdout.with_name('stderr.log')
    text = stdout.read_text() if stdout.exists() else ''
    progress = re.findall(r'Coordinate hashes verified (\d+) of (\d+)', text)
    counts = re.findall(r'Frozen log records (\d+) eligible sequence candidates (\d+)', text)
    result = dict(status='observed_original_full_afdb_catalog_and_readback_controllers',
        checked_utc=datetime.now(timezone.utc).isoformat(), original_handles=handles,
        source_plan_sha256=sha(pp), expected_taxa=526, expected_proteins=5815847,
        output_root_present=root.exists(), raw_catalog_receipt_present=(root/'receipt.json').exists(),
        producer_adapter_receipt_present=Path('metadata/completed_afdb_catalog_20261005_v1.json').exists(),
        independent_readback_receipt_present=Path('metadata/completed_afdb_catalog_readback_20261005_v1.json').exists(),
        last_logged_coordinate_hash_count=int(progress[-1][0]) if progress else None,
        logged_selected_models=int(progress[-1][1]) if progress else None,
        logged_inventory_records=int(counts[-1][0]) if counts else None,
        logged_eligible_sequence_candidates=int(counts[-1][1]) if counts else None,
        stderr_bytes_at_observation=stderr.stat().st_size if stderr.exists() else None,
        source_hashes={str(path):sha(path) for path in [pp, Path(__file__), *map(Path, launches)]},
        scientific_eligibility=False, new_downloads=0, new_predictions=0, gpu=False,
        scope='Exact original PID/create/command/cgroup or invocation-linked terminal observation. '
              'Logged candidate/model/hash counts and partial file presence are progress only; '
              'they do not prove completed coverage, every coordinate hash, independent ranking '
              'or original execution transport closure. All original taxa/proteins remain required.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['original_handles','source_hashes']}, indent=2))


if __name__ == '__main__':
    main()
