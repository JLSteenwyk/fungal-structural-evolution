#!/usr/bin/env python3
"""Recheck the complete historical census with direct mapped scalar lookup.

This reader uses a separate scalar traversal. Native reference calculations
and their double-roundtrip evidence are inherited; SciPy reference evaluation
shares the producer's declared function and is not a third numerical method.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

from ancestral_chain_attempt import sha
from assess_full_historical_baliphy_numbers_v2 import reference_gamma
from full_weighted_fit_exports import atomic
from independent_short_sampler_outputs import strict_json
from reference_measurement_union_sources import bind, verify


def scalar_value(row, key):
    if key == 'iter':return row['iter']
    if key in row['statistics//']:return row['statistics//'][key]
    value = row['parameters//'];parts = key.split('/')
    for group in parts[:-1]:value = value[group+'/']
    item = parts[-1]
    match = re.fullmatch(r'(.+)\[([^\[\]]+)\]', item)
    return value[match[1]][match[2]] if match else value[item]


def leaves(value):
    if isinstance(value, dict):return sum(leaves(x) for x in value.values())
    assert type(value) in [int, float] and math.isfinite(value)
    return 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--producer', type=Path, required=True);p.add_argument('--transport', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True);p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args();assert not a.receipt.exists()
    proof = strict_json(a.producer.read_text());transport = strict_json(a.transport.read_text())
    assert proof['status'] == 'completed_full_historical_baliphy_scalar_and_rate_encoding_assessment_v2'
    assert transport['validation_sha256'] == sha(a.producer) and transport['actual_tool_terminal_exit_code'] == 0
    verify(transport['source_hashes']);bindings = dict(transport['source_hashes'])
    for q in [a.producer, a.transport, Path(__file__)]:bind(bindings, q)
    assert proof['scientific_eligibility'] is proof['posterior_qualified'] is proof['original_numbers_repaired'] is False
    assert proof['native_mcmc_runs'] == 0 and proof['all_historical_rate_cells_unqualified'] == 388800
    native_plan = strict_json(Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json').read_text())
    dispositions = strict_json((Path(native_plan['output'])/'dispositions.json').read_text())
    assert len(dispositions) == len({r['chain_id'] for r in dispositions}) == 1620
    scalar_export = strict_json((a.output_root/'scalar-encoding-discrepancies.json').read_text())
    rate_export = strict_json((a.output_root/'rate-reference-assessment.json').read_text())
    frames = {r['id']: r for r in rate_export['frames']};assert len(frames) == len(rate_export['frames']) == 4860
    assert rate_export['native_alpha_parse_differences'] == [] and proof['native_alpha_parse_difference_count'] == 0
    native_path = Path(proof['native_receipt']);native_receipt = strict_json(native_path.read_text());assert native_receipt['exit_code'] == 0
    native_lines = (native_path.parent/'stdout.log').read_text().splitlines()
    json_lines = [];footer = []
    for line in native_lines:
        if line.startswith('{') and not footer:json_lines.append(line)
        else:footer.append(line)
    assert len(json_lines) == 4860 and len(footer) == 6 and footer[:2] == ['', 'Work:']
    assert footer[2].startswith('  start: ') and footer[3].startswith('    end: ')
    assert footer[4].startswith('  total (elapsed) time: ') and footer[5].startswith('  total (CPU) time: ')
    native_rows = {r['id']: r for r in map(strict_json, json_lines)};assert set(native_rows) == set(frames)
    counts = Counter();columns = Counter();families = Counter();priors = Counter();rate_counts = Counter();rate_families = Counter()
    discrepancies = [];invalid = [];native_roundtrips = 0;roles = 0
    for role in sorted(dispositions, key=lambda r: r['chain_id']):
        cid = role['chain_id'];assert role['scientific_eligibility'] is role['posterior_qualified'] is False
        assert role['status'] == 'full_short_sampler_output_integrity_checked_not_posterior'
        directory = next(Path(role['native_receipt']).parent.glob('independent-chain-*'))
        with (directory/'C1.log').open() as f:scalar_rows = list(csv.DictReader(f, delimiter='\t'))
        mapping = strict_json((directory/'C1.log.column-map.json').read_text())
        scalar_lines = (directory/'C1.log.json').read_text().splitlines()[1:]
        assert len(scalar_rows) == len(scalar_lines) == 21
        for iteration, (scalar, line) in enumerate(zip(scalar_rows, scalar_lines)):
            assert int(scalar['iter']) == iteration
            try:decoded = strict_json(line)
            except ValueError as error:
                assert line.startswith('{"iter":'+str(iteration)+',')
                invalid.append(dict(chain_id=cid, family=role['family'], prior=role['prior_label'],
                    iteration_from_mapped_tsv_and_original_json_prefix=iteration,
                    original_json_line_sha256=hashlib.sha256(line.encode()).hexdigest(),
                    source_json_path=str(directory/'C1.log.json'), error_type=type(error).__name__, error_message=str(error),
                    unavailable_numeric_comparisons=len(mapping), original_numbers_repaired=False, scientific_eligibility=False))
                counts['invalid_json_scalar_records_retained'] += 1;counts['unavailable_numeric_comparisons_invalid_json'] += len(mapping)
                continue
            assert set(decoded) == {'iter', 'statistics//', 'parameters//'} and decoded['iter'] == iteration
            assert leaves(decoded) == len(mapping)
            counts['valid_json_scalar_records_compared'] += 1
            for key, column in mapping.items():
                native_value, encoded = float(scalar[column]), scalar_value(decoded, key)
                assert math.isfinite(native_value)
                agrees = math.isclose(native_value, encoded, rel_tol=proof['scalar_relative_tolerance'], abs_tol=0.)
                counts['numeric_values_compared'] += 1
                counts['agrees_at_declared_relative_precision' if agrees else 'encoding_discrepancy_values'] += 1
                if not agrees:
                    discrepancies.append(dict(chain_id=cid, family=role['family'], prior=role['prior_label'],
                        iteration=iteration, original_json_key=key, mapped_tsv_column=column,
                        original_tsv_token=scalar[column], tsv_value=native_value, original_json_value=encoded,
                        relative_difference=abs(encoded-native_value)/abs(native_value) if native_value != 0 else None,
                        scientific_eligibility=False, original_numbers_repaired=False))
                    columns[column] += 1;families[role['family']] += 1;priors[role['prior_label']] += 1
        property_path = directory/'C1.P1.site-property-samples.jsonl'
        property_lines = property_path.read_text().splitlines();assert len(property_lines) == 3
        for iteration, line in zip([0, 10, 20], property_lines):
            frame = strict_json(line);assert frame['iter'] == iteration
            ident = cid+':'+str(iteration);record, native = frames[ident], native_rows[ident]
            assert record['chain_id'] == cid and record['iteration'] == iteration and record['family'] == role['family'] and record['prior'] == role['prior_label']
            rates = np.asarray(frame['properties']['rate']);assert rates.shape == (4, 20)
            assert np.array_equal(rates, np.repeat(rates[:, :1], 20, axis=1))
            original = rates[:, 0];actual = np.asarray(native['model_rates'])
            assert record['original_category_rates'] == original.tolist()
            assert record['original_property_path'] == str(property_path) and record['original_property_sha256'] == sha(property_path)
            token = scalar_rows[iteration]['ASRV.Gamma:alpha'];assert native['alpha_token'] == record['original_alpha_tsv_token'] == token
            alpha = float(token)
            assert record['alpha_tsv'] == (alpha if math.isfinite(alpha) else 'positive_infinity')
            assert record['pure_native_alpha'] == (native['alpha'] if math.isfinite(alpha) else 'positive_infinity')
            if math.isfinite(alpha):assert native['alpha'] == alpha
            else:assert alpha == math.inf and native['alpha'] is None
            assert record['pure_native_model_rates'] == actual.tolist() and native['roundtrips'] == [True]*4
            assert record['native_original_encoder_values'] == strict_json(native['original_encoding'])
            agreements = [math.isclose(float(x), float(y), rel_tol=proof['rate_relative_tolerance'], abs_tol=0.) for x, y in zip(original, actual)]
            old_agreements = [math.isclose(float(x), float(y), rel_tol=proof['rate_relative_tolerance'], abs_tol=0.) for x, y in zip(original, record['native_original_encoder_values'])]
            assert record['original_vs_native_category_agreement'] == agreements
            assert record['original_vs_native_original_encoder_agreement'] == old_agreements
            mean_pass = math.isclose(math.fsum(original)/4, 1., rel_tol=1e-10, abs_tol=1e-10)
            assert record['mean_one_pass'] is mean_pass
            assert record['historical_rate_properties_qualified'] is record['original_numbers_repaired'] is record['scientific_eligibility'] is False
            independent, route = reference_gamma(alpha)
            check = dict(rates=independent.tolist(), route=route,
                native_agrees_at_declared_precision=bool(np.allclose(independent, actual, rtol=proof['independent_gamma_relative_tolerance'], atol=proof['independent_gamma_absolute_tolerance'])),
                maximum_absolute_difference=float(np.max(abs(independent-actual))))
            assert record['independent_gamma_reference'] == check and record['independent_reference_qualified'] is False
            changed = 4-sum(agreements)
            rate_counts['category_values_compared'] += 4;rate_counts['repeated_amino_acid_rate_cells_compared'] += 80
            rate_counts['rate_discrepancy_categories'] += changed;rate_counts['rate_discrepancy_repeated_cells'] += 20*changed
            rate_counts['mean_one_pass_frames' if mean_pass else 'mean_one_failed_frames'] += 1
            if changed:
                rate_counts['rate_discrepancy_frames'] += 1;rate_families[role['family']] += 1
                if mean_pass:rate_counts['rate_discrepancy_frames_passing_mean_one'] += 1
                if all(old_agreements):rate_counts['discrepancy_frames_consistent_with_original_encoder'] += 1
            rate_counts['independent_gamma_agreement_frames' if check['native_agrees_at_declared_precision'] else 'independent_gamma_precision_review_frames'] += 1
            if not math.isfinite(alpha):rate_counts['positive_infinite_alpha_reference_frames'] += 1
            native_roundtrips += 4
        roles += 1
        if roles % 100 == 0:print('historical_numeric_independent_lookup_readback', roles, '/1620', flush=True)
    assert scalar_export['rows'] == discrepancies and scalar_export['invalid_scalar_records'] == invalid
    assert scalar_export['scalar_summary'] == proof['scalar_numeric_counts'] == dict(counts)
    assert proof['scalar_discrepancies_by_column'] == dict(columns) and proof['scalar_discrepancies_by_family'] == dict(families) and proof['scalar_discrepancies_by_prior'] == dict(priors)
    assert rate_export['counts'] == proof['rate_counts'] == dict(rate_counts)
    assert proof['rate_discrepancy_frames_by_family'] == dict(rate_families)
    assert roles == 1620 and native_roundtrips == 19440
    verify(bindings)
    result = dict(status='passed_full_historical_number_assessment_direct_lookup_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), producer_receipt_sha256=sha(a.producer),
        full_roles=1620, full_property_frames=4860, scalar_rows_compared=34020,
        scalar_numeric_counts=dict(counts), rate_counts=dict(rate_counts),
        native_double_roundtrips=native_roundtrips, native_roundtrip_checks_inherited=True,
        separate_direct_scalar_lookup=True, native_rate_probe_repeated=False,
        scipy_reference_function_shared=True, original_numbers_repaired=False, native_mcmc_runs=0,
        all_historical_rate_cells_unqualified=388800, scientific_eligibility=False, posterior_qualified=False,
        gpu=False, source_hashes=bindings,
        scope='All original mapped TSV/JSON rows and property values reconstructed with a separate direct scalar lookup; '
            'every discrepancy/malformed-record/status/count and full serialized rate/alpha/reference field checked. '
            'Native batch/roundtrip and original-source closure evidence inherited; SciPy function is shared. '
            'No restored native historical floats, repaired output, posterior acceptance, new MCMC or installed correction.')
    atomic(a.receipt, result)
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__':main()
