#!/usr/bin/env python3
"""Check full-reader integration on closed native traces and incomplete/error states."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from readback_baliphy_log_alpha_v10_full_grid_v1 import diagnostics
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    args.output.mkdir(exist_ok=False)
    closed_path = Path('metadata/baliphy_log_alpha_v10_completed_20261005_v1.json')
    qualification_path = Path('metadata/baliphy_log_alpha_v10_software_validation_20261005_v1.json')
    closed, qualified = [json.loads(p.read_text()) for p in [closed_path, qualification_path]]
    assert closed['status'] == 'complete_verified_V10_latent_alpha_logger_software_and_diagnostic_replay'
    verify(closed['source_hashes'])
    pins = dict(closed['source_hashes'])
    positives = []
    for outcome in qualified['outcomes']:
        current = Path(outcome['native_receipt']).parent / 'independent-chain-1'
        result = diagnostics(current / 'C1.P1.log-alpha-samples.jsonl', current / 'C1.log.json', outcome['prior'], True)
        assert result['checked_rows'] == result['diagnostic_rows'] == 21
        assert result['all_saved_rows_checked'] and result['complete_native_trace_checked']
        assert not result['errors'] and result['overflow_rows'] == 0
        positives.append(dict(prior=outcome['prior'], **result))
    assert len(positives) == 3 and sum(r['checked_rows'] for r in positives) == 63
    first = qualified['outcomes'][0]
    current = Path(first['native_receipt']).parent / 'independent-chain-1'
    diagnostic = (current / 'C1.P1.log-alpha-samples.jsonl').read_text().splitlines()
    reference = (current / 'C1.log.json').read_text().splitlines()
    altered = list(diagnostic)
    value = json.loads(altered[1])
    value['parameters//']['S1/']['latentLogDensity'] += 1
    altered[1] = json.dumps(value)
    bad = args.output / 'changed_density.jsonl'
    bad.write_text('\n'.join(altered) + '\n')
    result = diagnostics(bad, current / 'C1.log.json', first['prior'], True)
    assert result['diagnostic_rows'] == 21 and result['checked_rows'] == 20 and len(result['errors']) == 1
    assert not result['complete_native_trace_checked']
    short = args.output / 'partial_diagnostic.jsonl'
    short_ref = args.output / 'partial_reference.jsonl'
    short.write_text('\n'.join(diagnostic[:-1]) + '\n')
    short_ref.write_text('\n'.join(reference[:-1]) + '\n')
    result = diagnostics(short, short_ref, first['prior'], True)
    assert result['checked_rows'] == 20 and result['all_saved_rows_checked'] and not result['complete_native_trace_checked']
    result = diagnostics(current / 'C1.P1.log-alpha-samples.jsonl', current / 'C1.log.json', first['prior'], False)
    assert result['checked_rows'] == 21 and result['all_saved_rows_checked'] and not result['complete_native_trace_checked']
    result = diagnostics(args.output / 'missing.jsonl', current / 'C1.log.json', first['prior'], True)
    assert result['checked_rows'] == 0 and result['errors'] == [dict(kind='missing_trace')]
    controls = ['changed_density_retained_as_review', 'partial_trace_not_complete',
                'failed_native_not_complete_despite_full_trace', 'missing_diagnostic_retained_as_review']
    for path in [Path(__file__), closed_path, qualification_path,
                 *[Path('scripts', name + '.py') for name in ['readback_baliphy_log_alpha_v10_full_grid_v1',
                    'readback_baliphy_log_alpha_v10', 'run_baliphy_log_alpha_v10_full_grid_v1',
                    'baliphy_log_alpha_v10_full_pairs_v1', 'baliphy_log_alpha_v10_full_jobs_v1']],
                 bad, short, short_ref]:
        bind(pins, path)
    verify(pins)
    result = dict(status='passed_V10_full_reader_closed_trace_and_review_state_integration',
        checked_utc=datetime.now(timezone.utc).isoformat(), native_controls_reused=3,
        saved_native_rows=63, review_controls=controls, positive_results=positives,
        source_hashes=pins, scientific_eligibility=False, posterior_qualified=False,
        new_native_runs=0, gpu=False,
        scope='Integration checks use all63closed native V10control rows and four deliberately incomplete/error states. '
              'No new native model or corpus subset; complete24role full-input run/readback remains required.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'positive_results']}, indent=2))


if __name__ == '__main__':
    main()
