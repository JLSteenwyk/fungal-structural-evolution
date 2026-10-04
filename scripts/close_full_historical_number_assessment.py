#!/usr/bin/env python3
"""Bind full historical assessment/readback to both original driver journals."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from full_weighted_fit_exports import atomic
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True);a = p.parse_args()
    paths = ['metadata/full_historical_baliphy_number_assessment_20261004_v2.json',
        'metadata/full_historical_baliphy_number_readback_20261004_v1.json']
    transports = ['metadata/full_historical_baliphy_number_assessment_transport_20261004_v2.json',
        'metadata/full_historical_baliphy_number_readback_transport_20261004_v1.json']
    producer, reader = [json.loads(Path(n).read_text()) for n in paths]
    assert producer['status'] == 'completed_full_historical_baliphy_scalar_and_rate_encoding_assessment_v2'
    assert reader['status'] == 'passed_full_historical_number_assessment_direct_lookup_readback'
    assert reader['producer_receipt_sha256'] == sha(paths[0])
    for record in [producer, reader]:
        assert (record['full_roles'], record['full_property_frames'], record['scalar_rows_compared'], record['native_double_roundtrips']) == (1620, 4860, 34020, 19440)
        assert record['scientific_eligibility'] is record['posterior_qualified'] is record['original_numbers_repaired'] is False
        assert record['all_historical_rate_cells_unqualified'] == 388800 and record['native_mcmc_runs'] == 0
    assert producer['scalar_numeric_counts'] == reader['scalar_numeric_counts']
    assert producer['rate_counts'] == reader['rate_counts']
    assert reader['separate_direct_scalar_lookup'] is True and reader['native_roundtrip_checks_inherited'] is True
    bindings = {};drivers = []
    for receipt_path, transport_path in zip(paths, transports):
        t = json.loads(Path(transport_path).read_text());verify(t['source_hashes'])
        assert t['actual_tool_terminal_exit_code'] == 0 and t['validation_sha256'] == sha(receipt_path)
        assert t['entire_terminal_payload_matched'] is True
        assert t['original_start_records'] == t['original_completion_records'] == 1
        assert t['exact_wrapper_pid_journal_entries'] == 2
        for name, d in t['source_hashes'].items():bind(bindings, name, d)
        bind(bindings, transport_path);bind(bindings, receipt_path)
        drivers.append(dict(unit=t['unit'], invocation_id=t['invocation_id'], wrapper=t['wrapper'],
            actual_tool_session_id=t['actual_tool_session_id'], actual_tool_terminal_exit_code=0,
            transport=transport_path, transport_sha256=sha(transport_path)))
    assert len({d['invocation_id'] for d in drivers}) == 2
    bind(bindings, Path(__file__));verify(bindings)
    result = dict(status='complete_verified_full_historical_scalar_and_rate_encoding_assessment',
        checked_utc=datetime.now(timezone.utc).isoformat(), full_roles=1620, full_model_groups=405,
        full_effective_inputs=135, full_configuration_aliases=324, scalar_rows_assessed=34020,
        full_property_frames=4860, scalar_numeric_counts=producer['scalar_numeric_counts'],
        rate_counts=producer['rate_counts'], all_historical_rate_cells_unqualified=388800,
        producer_receipt=paths[0], producer_receipt_sha256=sha(paths[0]),
        independent_readback=paths[1], independent_readback_sha256=sha(paths[1]),
        original_drivers=drivers, exact_original_driver_journals_checked=2,
        source_hashes=bindings, bound_source_hashes=len(bindings),
        scientific_eligibility=False, posterior_qualified=False, original_numbers_repaired=False,
        installed_formatter_fixed=False, native_mcmc_runs=0, gpu=False, new_cost_usd=0,
        scope='Full historical scalar/rate forensic accounting closes with separate direct scalar lookup, '
            'all original outcomes and both exact original driver waits/journal terminal payloads. '
            'Original source closure and native reference/roundtrip facts inherited, all consumed bytes hashed. '
            'Not a complete census of unobservable original float loss, recovery of exact historical native values, '
            'installed formatter/crash repair, posterior/biological acceptance or completion of any evolutionary aim.')
    atomic(a.output, result)
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'original_drivers']}, indent=2))


if __name__ == '__main__':main()
