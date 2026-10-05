#!/usr/bin/env python3
"""Replay all V10 paired files and latent states before serialized native readback."""
import argparse
from datetime import datetime, timezone
from decimal import Decimal
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from readback_baliphy_log_alpha_v10 import load, check, HEADER
from reference_measurement_union_sources import bind, verify
from run_baliphy_log_alpha_v10_full_grid_v1 import run


FILES = ['C1.log', 'C1.log.json', 'C1.log.column-map.json', 'runtime-tree.nwk',
         'C1.P1.fastas', 'C1.P1.site-property-samples.jsonl']


def equal_bytes(first, second):
    with first.open('rb') as left, second.open('rb') as right:
        while True:
            a = left.read(32768)
            b = right.read(32768)
            if a != b:
                return False
            if not a:
                return True


def diagnostics(path, reference_path, prior, native_zero):
    result = dict(diagnostic_rows=0, checked_rows=0, overflow_rows=0, errors=[],
                  overflow_observations=[], all_saved_rows_checked=False,
                  complete_native_trace_checked=False, scientific_eligibility=False)
    if not path.is_file() or not reference_path.is_file():
        result['errors'].append(dict(kind='missing_trace'))
        return result
    try:
        rows = [load(line) for line in path.read_text().splitlines()]
        original = [load(line) for line in reference_path.read_text().splitlines()]
        assert rows[0] == original[0] == HEADER
        rows, original = rows[1:], original[1:]
        result['diagnostic_rows'] = len(rows)
        assert len(rows) == len(original)
        assert [r['iter'] for r in rows] == list(range(len(rows)))
    except (AssertionError, ValueError, KeyError, TypeError, OSError) as error:
        result['errors'].append(dict(kind='trace_structure_review', error_type=type(error).__name__, error=str(error)))
        return result
    mu, scale = {'broad': (Decimal(0), Decimal(2)), 'centered': (Decimal(0), Decimal(1)),
                 'package': (Decimal(6), Decimal(2))}[prior]
    for row, reference in zip(rows, original):
        try:
            overflow = check(row, mu, scale, reference)
            result['checked_rows'] += 1
            result['overflow_rows'] += int(overflow)
            if overflow:
                state = row['parameters//']['S1/']
                result['overflow_observations'].append(dict(iteration=row['iter'],
                    finite_latent_log_alpha=str(state['latentLogAlpha']),
                    latent_log_density=str(state['latentLogDensity']),
                    derived_alpha=state['derivedAlpha'], category_rates=[str(v) for v in state['categoryRates']]))
        except (AssertionError, ValueError, KeyError, TypeError, ArithmeticError) as error:
            result['errors'].append(dict(kind='latent_numeric_review', iteration=row.get('iter'),
                error_type=type(error).__name__, error=str(error)))
    result['all_saved_rows_checked'] = not result['errors'] and result['checked_rows'] == len(rows)
    result['complete_native_trace_checked'] = native_zero and len(rows) == 21 and result['all_saved_rows_checked']
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['plan', 'producer-receipt', 'producer-transport', 'receipt']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    plan, producer, transport = [json.loads(p.read_text()) for p in
                                [args.plan, args.producer_receipt, args.producer_transport]]
    assert producer['status'] == 'complete_all24_V10_full_comparison_dispositions_pending_readback'
    assert transport['original_tool_terminal_exit_code'] == 0
    assert transport['validation_sha256'] == sha(args.producer_receipt)
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    verify(transport['source_hashes'])
    root = Path(plan['output'])
    jobs = json.loads(Path(plan['jobs']).read_text())
    rows = {r['chain_id']: r for r in json.loads((root / 'dispositions.json').read_text())}
    saved_pairs = json.loads((root / 'paired_scientific_outputs.json').read_text())
    assert len(jobs) == len(rows) == len(saved_pairs) == 24
    pairs, audits = [], []
    for job in jobs:
        cid = job['chain']['chain_id']
        row = rows[cid]
        old_receipt = Path(job['paired_v7_native_receipt'])
        assert sha(old_receipt) == job['paired_v7_native_receipt_sha256']
        old, current = [p / 'independent-chain-1' for p in
                        [old_receipt.parent, Path(row['native_receipt']).parent]]
        files = []
        for name in FILES:
            first, second = old / name, current / name
            assert first.is_file()
            files.append(dict(name=name, original_path=str(first), current_path=str(second),
                original_sha256=sha(first), current_sha256=sha(second) if second.is_file() else None,
                current_present=second.is_file(), byte_identical=second.is_file() and equal_bytes(first, second)))
        path = current / 'C1.P1.log-alpha-samples.jsonl'
        count = max(0, len(path.read_text().splitlines()) - 1) if path.is_file() else 0
        identical = row['exit_code'] == 0 and all(f['byte_identical'] for f in files)
        pairs.append(dict(chain_id=cid, source_v7_chain_id=job['source_v7_chain_id'], native_exit_code=row['exit_code'],
            files=files, all_original_files_byte_identical=identical, diagnostic_path=str(path),
            diagnostic_present=path.is_file(), diagnostic_sha256=sha(path) if path.is_file() else None,
            diagnostic_rows_present=count,
            status='all_original_scientific_files_byte_identical_pending_diagnostic_readback' if identical else
                   'native_or_paired_scientific_output_review_required', scientific_eligibility=False, posterior_qualified=False))
        audits.append(dict(chain_id=cid, **diagnostics(path, current / 'C1.log.json', job['chain']['prior_label'], row['exit_code'] == 0)))
    assert sorted(pairs, key=lambda r: r['chain_id']) == saved_pairs
    assert producer['paired_roles_with_identical_scientific_files'] == sum(r['all_original_files_byte_identical'] for r in pairs)
    native = run(args.plan, reader=True)
    pins = dict(native['source_hashes'])
    for path, digest in transport['source_hashes'].items():
        bind(pins, path, digest)
    for path in [args.plan, args.producer_receipt, args.producer_transport, Path(__file__),
                 Path('scripts/readback_baliphy_log_alpha_v10.py'), root / 'paired_scientific_outputs.json', root / 'readback.json']:
        bind(pins, path)
    verify(pins)
    result = dict(status='complete_all24_V10_independent_paired_files_latent_states_and_native_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), roles=24, paired_file_comparisons=sum(len(r['files']) for r in pairs),
        paired_roles_with_identical_scientific_files=sum(r['all_original_files_byte_identical'] for r in pairs),
        diagnostic_rows=sum(r['diagnostic_rows'] for r in audits), checked_rows=sum(r['checked_rows'] for r in audits),
        overflow_rows=sum(r['overflow_rows'] for r in audits), diagnostic_audits=audits,
        complete_diagnostic_roles=sum(r['complete_native_trace_checked'] for r in audits),
        full_input_logger_noninterference_passed=all(r['all_original_files_byte_identical'] for r in pairs),
        full_latent_trace_integrity_passed=all(r['complete_native_trace_checked'] for r in audits),
        native_summary={k: native[k] for k in ['native_zero_exit_roles', 'integrity_checked_roles', 'saved_joint_frames']},
        source_hashes=pins, scientific_eligibility=False, posterior_qualified=False, gpu=False,
        scope='All24paired native full-input roles retained; separate byte comparison and qualified strict JSON/90digit Decimal '
              'reader check every available latent row, context, derived alpha, density/rates/quality state. Errors and native '
              'failures remain review outcomes; no automatic retry, tolerance relaxation, old-array admission or adequate posterior '
              'claim. Existing native/scalar/joint serialized inspectors and job contracts are shared dependencies. New captures '
              'can explain these newly observed states, not recover missing historical latent values.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'diagnostic_audits']}, indent=2))


if __name__ == '__main__':
    main()
