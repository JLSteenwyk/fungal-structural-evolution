#!/usr/bin/env python3
"""Read-only scalar/rate encoding assessment of every historical short role.

No damaged value is repaired, and no passing comparison qualifies a historical
number. Pure native rate calculations use logged TSV alpha; their outputs are
diagnostic references, not recovered original binary64 sampler state.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammainc, gammaincc, gammaincinv
from scipy.stats import norm

from ancestral_chain_attempt import run_attempt, sha
from full_exact_covariance_sources import closed_subset
from full_weighted_fit_exports import atomic
from independent_short_sampler_outputs import strict_json
from reference_measurement_union_sources import bind, verify

SCALAR_RTOL = RATE_RTOL = 2e-13
NATIVE_SCIPY_RTOL = 1e-7
NATIVE_SCIPY_ATOL = 2e-12


def flatten_scalar(row):
    if set(row) != {'iter', 'statistics//', 'parameters//'}:
        raise ValueError('Unexpected scalar record schema')
    result = {'iter': row['iter']}
    def add(name, value):
        if name in result or type(value) not in [int, float] or not math.isfinite(value):
            raise ValueError('Duplicate/nonfinite/nonnumeric scalar field')
        result[name] = value
    def walk(fields, prefix):
        for key, value in fields.items():
            if not isinstance(key, str):raise ValueError('Scalar field name required')
            if isinstance(value, dict):
                if key.endswith('/'):
                    walk(value, prefix+key)
                else:
                    for subkey, item in value.items():add(prefix+key+'['+subkey+']', item)
            else:add(prefix+key, value)
    walk(row['statistics//'], '');walk(row['parameters//'], '')
    return result


def close(a, b, rtol):
    # Zero values are exact; no absolute tolerance hides tiny exponent damage.
    return math.isclose(a, b, rel_tol=rtol, abs_tol=0.)


def reference_gamma(alpha):
    if alpha == math.inf:return np.ones(4), 'positive_infinite_alpha_equal_rate_limit_reference'
    if not math.isfinite(alpha) or alpha <= 0:raise ValueError('Finite positive TSV alpha required')
    if alpha+1 == alpha:
        densities = np.r_[0., norm.pdf(norm.ppf([.25, .5, .75])), 0.]
        rates = 1.+4./np.sqrt(alpha)*np.diff(-densities)
        route = 'large_alpha_normal_limit_reference'
    else:
        boundaries = gammaincinv(alpha, [.25, .5, .75])
        lower = np.r_[0., gammainc(alpha+1., boundaries), 1.]
        upper = np.r_[1., gammaincc(alpha+1., boundaries), 0.]
        rates = np.r_[4*np.diff(lower)[:2], -4*np.diff(upper)[2:]]
        route = 'gamma_conditional_bin_mean_reference'
    if not np.isfinite(rates).all() or np.any(rates < 0) or not rates.sum() > 0:
        raise ArithmeticError('Independent gamma reference requires precision review')
    return rates/(rates.sum()/4), route


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True);p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args();assert not a.receipt.exists();a.output.mkdir(parents=True, exist_ok=False)
    bindings = {};native_plan = Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json')
    plan = strict_json(native_plan.read_text());root = Path(plan['output'])
    dispositions_path = root/'dispositions.json';jobs_path = Path(plan['jobs'])
    consumed = [native_plan, dispositions_path, jobs_path, root/'receipt.json', root/'readback.json']
    completed = closed_subset(plan['completion'], 'complete_verified_full_reference_short_sampler_qualification', consumed, bindings)
    assert completed['full_chains'] == completed['checked_sampler_attempts'] == 1620
    assert completed['full_quartets'] == 405 and completed['saved_alignments'] == 4860
    archive = strict_json(Path(completed['full_hash_archive']).read_text())
    original_hashes = archive['source_hashes']
    jobs = strict_json(jobs_path.read_text());dispositions = strict_json(dispositions_path.read_text())
    job_by_id = {j['chain']['chain_id']: j for j in jobs};role_by_id = {r['chain_id']: r for r in dispositions}
    assert len(jobs) == len(dispositions) == len(job_by_id) == len(role_by_id) == 1620
    assert set(job_by_id) == set(role_by_id)
    for name in ['assess_full_historical_baliphy_numbers.py', 'ancestral_chain_attempt.py', 'full_exact_covariance_sources.py']:
        bind(bindings, 'scripts/'+name)
    probe_source = Path('config/native-format-probes/full-gamma-rate-reference-v1.hs');bind(bindings, probe_source)
    binary = Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/bali-phy').resolve()
    api = binary.parent.parent/'lib/bali-phy/haskell'
    libraries = [api/name for name in ['SModel/ASRV.hs', 'SModel/MixtureModel.hs',
        'Probability/Distribution/Discrete.hs', 'Probability/Logger.hs', 'Data/JSON.hs',
        'Data/JSON/Encoding.hs', 'Data/JSON/Types/ToJSON.hs', 'Data/JSON/Types/Foreign.hs', 'Data/Text/IO.hs']]
    for q in libraries:bind(bindings, q)
    primary_evidence = Path('metadata/baliphy_native_exponent_formatting_failure_20261003_v1.json')
    primary = strict_json(primary_evidence.read_text());bind(bindings, primary_evidence)
    for source in primary['primary_sources']:bind(bindings, source['path'], source['sha256'])
    scalar_discrepancies = [];invalid_scalar_records = [];scalar_summary = Counter();column_counts = Counter();family_scalar = Counter();prior_scalar = Counter()
    frames = [];input_lines = [];frame_role = {};row_counts = Counter()
    for cid in sorted(job_by_id):
        job, role = job_by_id[cid], role_by_id[cid]
        assert role['status'] == 'full_short_sampler_output_integrity_checked_not_posterior'
        assert role['exit_code'] == 0 and role['seed'] == job['chain']['seed']
        receipt_path = Path(role['native_receipt']);bind(bindings, receipt_path, role['native_receipt_sha256'])
        receipt = strict_json(receipt_path.read_text());assert receipt['exit_code'] == 0
        directories = list(receipt_path.parent.glob('independent-chain-*'));assert len(directories) == 1
        directory = directories[0]
        paths = [directory/name for name in ['C1.log', 'C1.log.json', 'C1.log.column-map.json', 'C1.P1.site-property-samples.jsonl']]
        for path in paths:
            assert str(path) in original_hashes
            assert receipt['artifacts'][str(path.relative_to(receipt_path.parent))] == original_hashes[str(path)]
            bind(bindings, path, original_hashes[str(path)])
        log, jsonlog, mapping_path, property_path = paths
        mapping = strict_json(mapping_path.read_text())
        assert len(set(mapping.values())) == len(mapping)
        with log.open() as f:
            reader = csv.DictReader(f, delimiter='\t');header = reader.fieldnames;rows = list(reader)
        assert len(header) == len(set(header)) and set(mapping.values()) == set(header)
        assert [int(r['iter']) for r in rows] == list(range(21))
        json_lines = jsonlog.read_text().splitlines();metadata = strict_json(json_lines[0])
        assert metadata == dict(fields=['iter', 'prior', 'likelihood', 'posterior'], nested=True, format='MCON', version='0.2')
        assert len(json_lines) == 22
        for scalar, line in zip(rows, json_lines[1:]):
            iteration = int(scalar['iter'])
            try:decoded = flatten_scalar(strict_json(line))
            except ValueError as error:
                assert line.startswith('{"iter":'+str(iteration)+',')
                invalid_scalar_records.append(dict(chain_id=cid, family=role['family'], prior=role['prior_label'],
                    iteration_from_mapped_tsv_and_original_json_prefix=iteration,
                    original_json_line_sha256=hashlib.sha256(line.encode()).hexdigest(),
                    source_json_path=str(jsonlog), error_type=type(error).__name__, error_message=str(error),
                    unavailable_numeric_comparisons=len(mapping), original_numbers_repaired=False, scientific_eligibility=False))
                scalar_summary['invalid_json_scalar_records_retained'] += 1
                scalar_summary['unavailable_numeric_comparisons_invalid_json'] += len(mapping)
                continue
            assert set(decoded) == set(mapping) and decoded['iter'] == iteration
            scalar_summary['valid_json_scalar_records_compared'] += 1
            for key, column in mapping.items():
                native_value = float(scalar[column]);json_value = decoded[key]
                assert math.isfinite(native_value)
                agrees = close(native_value, json_value, SCALAR_RTOL)
                scalar_summary['numeric_values_compared'] += 1
                scalar_summary['agrees_at_declared_relative_precision' if agrees else 'encoding_discrepancy_values'] += 1
                if not agrees:
                    scalar_discrepancies.append(dict(chain_id=cid, family=role['family'], prior=role['prior_label'],
                        iteration=iteration, original_json_key=key, mapped_tsv_column=column,
                        original_tsv_token=scalar[column], tsv_value=native_value, original_json_value=json_value,
                        relative_difference=abs(json_value-native_value)/abs(native_value) if native_value != 0 else None,
                        scientific_eligibility=False, original_numbers_repaired=False))
                    column_counts[column] += 1;family_scalar[role['family']] += 1;prior_scalar[role['prior_label']] += 1
        row_counts['scalar_rows'] += 21
        property_lines = property_path.read_text().splitlines();assert len(property_lines) == 3
        for iteration, line in zip([0, 10, 20], property_lines):
            frame = strict_json(line);assert frame['iter'] == iteration
            rates = np.asarray(frame['properties']['rate'], dtype=float)
            assert rates.shape == (4, 20) and np.isfinite(rates).all() and np.all(rates >= 0)
            assert np.array_equal(rates, np.repeat(rates[:, :1], 20, axis=1))
            token = rows[iteration]['ASRV.Gamma:alpha'];alpha = float(token)
            assert alpha > 0 and (math.isfinite(alpha) or alpha == math.inf)
            ident = cid+':'+str(iteration);assert ident not in frame_role
            input_lines.append(ident+' '+token+'\n')
            record = dict(id=ident, chain_id=cid, family=role['family'], prior=role['prior_label'], iteration=iteration,
                original_alpha_tsv_token=token, alpha_tsv=alpha if math.isfinite(alpha) else 'positive_infinity',
                original_category_rates=rates[:, 0].tolist(),
                original_property_path=str(property_path), original_property_sha256=sha(property_path),
                mean_one_pass=math.isclose(float(rates[:, 0].sum()/4), 1., rel_tol=1e-10, abs_tol=1e-10),
                historical_rate_properties_qualified=False, original_numbers_repaired=False, scientific_eligibility=False)
            frames.append(record);frame_role[ident] = record
        if len(frames) % 300 == 0:print('historical_numeric_source_frames', len(frames), '/4860', flush=True)
    assert len(frames) == 4860 and row_counts['scalar_rows'] == 34020
    inputs_path = a.output.resolve()/'logged-alphas.tsv'
    with inputs_path.open('x') as f:f.writelines(input_lines)
    native_program = a.output.resolve()/'rate-reference.hs'
    text = probe_source.read_text();assert text.count('__INPUT_PATH__') == 1
    native_program.write_text(text.replace('__INPUT_PATH__', str(inputs_path)))
    pin_paths = [Path('/usr/bin/prlimit'), binary, inputs_path, native_program, *libraries]
    config = dict(command=[str(pin_paths[0]), '--as='+str(8*2**30), '--cpu=600', '--fsize='+str(64*2**20),
        '--', str(binary), 'run', str(native_program)], timeout_seconds=900,
        pins={str(path): sha(path) for path in pin_paths})
    native_receipt = run_attempt(a.output/'native-reference', config)
    native = strict_json(native_receipt.read_text());assert native['exit_code'] == 0
    rows = [strict_json(line) for line in (native_receipt.parent/'stdout.log').read_text().splitlines()]
    assert len(rows) == len({r['id'] for r in rows}) == 4860 and {r['id'] for r in rows} == set(frame_role)
    for name, d in config['pins'].items():bind(bindings, name, d)
    for name, d in native['artifacts'].items():bind(bindings, native_receipt.parent/name, d)
    bind(bindings, native_receipt);bind(bindings, native_receipt.parent.parent/'configuration.json')
    rate_counts = Counter();rate_family = Counter();roundtrips = 0;alpha_read_differences = []
    for row in rows:
        record = frame_role[row['id']];original = np.asarray(record['original_category_rates'])
        assert row['alpha_token'] == record['original_alpha_tsv_token']
        if record['alpha_tsv'] == 'positive_infinity':
            assert row['alpha'] is None
            alpha = math.inf
            rate_counts['positive_infinite_alpha_reference_frames'] += 1
        else:
            alpha = row['alpha'];assert close(alpha, record['alpha_tsv'], 2e-15)
            if alpha != record['alpha_tsv']:alpha_read_differences.append(dict(id=row['id'], native_parsed_alpha=alpha, tsv_python_alpha=record['alpha_tsv']))
        rates = np.asarray(row['model_rates']);assert rates.shape == (4,) and np.isfinite(rates).all() and np.all(rates >= 0)
        assert row['roundtrips'] == [True]*4;roundtrips += 4
        assert math.isclose(float(rates.sum()/4), 1., rel_tol=1e-10, abs_tol=1e-10)
        agrees = [close(float(x), float(y), RATE_RTOL) for x, y in zip(original, rates)]
        old_encoding = strict_json(row['original_encoding']);assert len(old_encoding) == 4
        old_consistent = [close(float(x), float(y), RATE_RTOL) for x, y in zip(original, old_encoding)]
        try:
            reference, route = reference_gamma(alpha)
            independent_pass = bool(np.allclose(reference, rates, rtol=NATIVE_SCIPY_RTOL, atol=NATIVE_SCIPY_ATOL))
            reference_result = dict(rates=reference.tolist(), route=route,
                native_agrees_at_declared_precision=independent_pass,
                maximum_absolute_difference=float(np.max(abs(reference-rates))))
        except (ValueError, ArithmeticError) as error:
            reference_result = dict(error_type=type(error).__name__, error_message=str(error), native_agrees_at_declared_precision=False)
        record.update(pure_native_alpha=alpha if math.isfinite(alpha) else 'positive_infinity', pure_native_model_rates=rates.tolist(),
            native_original_encoder_values=old_encoding, original_vs_native_category_agreement=agrees,
            original_vs_native_original_encoder_agreement=old_consistent,
            independent_gamma_reference=reference_result, independent_reference_qualified=False)
        discrepancies = 4-sum(agrees)
        rate_counts['category_values_compared'] += 4;rate_counts['repeated_amino_acid_rate_cells_compared'] += 80
        rate_counts['rate_discrepancy_categories'] += discrepancies
        rate_counts['rate_discrepancy_repeated_cells'] += 20*discrepancies
        rate_counts['mean_one_pass_frames' if record['mean_one_pass'] else 'mean_one_failed_frames'] += 1
        if discrepancies:
            rate_counts['rate_discrepancy_frames'] += 1;rate_family[record['family']] += 1
            if record['mean_one_pass']:rate_counts['rate_discrepancy_frames_passing_mean_one'] += 1
            if all(old_consistent):rate_counts['discrepancy_frames_consistent_with_original_encoder'] += 1
        rate_counts['independent_gamma_agreement_frames' if reference_result['native_agrees_at_declared_precision'] else 'independent_gamma_precision_review_frames'] += 1
    assert rate_counts['mean_one_failed_frames'] == 12 and roundtrips == 19440
    assert scalar_summary['valid_json_scalar_records_compared']+scalar_summary['invalid_json_scalar_records_retained'] == 34020
    atomic(a.output/'scalar-encoding-discrepancies.json', dict(rows=scalar_discrepancies,
        invalid_scalar_records=invalid_scalar_records, scalar_summary=dict(scalar_summary)))
    atomic(a.output/'rate-reference-assessment.json', dict(frames=frames, counts=dict(rate_counts), native_alpha_parse_differences=alpha_read_differences))
    for q in a.output.rglob('*'):
        if q.is_file():bind(bindings, q)
    verify(bindings)
    result = dict(status='completed_full_historical_baliphy_scalar_and_rate_encoding_assessment_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(), full_roles=1620, full_model_groups=405,
        full_effective_inputs=135, full_configuration_aliases=324, scalar_rows_compared=34020,
        scalar_numeric_counts=dict(scalar_summary), scalar_discrepancies_by_column=dict(column_counts),
        scalar_discrepancies_by_family=dict(family_scalar), scalar_discrepancies_by_prior=dict(prior_scalar),
        full_property_frames=4860, rate_counts=dict(rate_counts), rate_discrepancy_frames_by_family=dict(rate_family),
        native_double_roundtrips=19440, native_alpha_parse_difference_count=len(alpha_read_differences),
        native_receipt=str(native_receipt), scalar_relative_tolerance=SCALAR_RTOL, scalar_absolute_tolerance=0.,
        rate_relative_tolerance=RATE_RTOL, rate_absolute_tolerance=0.,
        independent_gamma_relative_tolerance=NATIVE_SCIPY_RTOL, independent_gamma_absolute_tolerance=NATIVE_SCIPY_ATOL,
        all_historical_rate_cells_unqualified=388800, original_numbers_repaired=False,
        original_sources_or_outputs_changed=False, native_mcmc_runs=0, posterior_qualified=False,
        scientific_eligibility=False, gpu=False, new_cost_usd=0, source_hashes=bindings,
        scope='Full historical1620-role scalar JSON/TSV mapped-number and4860property-frame gamma reference assessment. '
            'Actual pure installed native rates at logged TSV alpha, CJSON native roundtrip and original-encoder '
            'diagnostics supplemented by SciPy conditional-bin-mean/large-alpha references; numeric reviews retained. '
            'This does not recover exact original native floats, qualify historical rates, repair scalar JSON, '
            'prove likelihood/mixing or fix installed formatter/crashes. Every original artifact remains unchanged; '
            'original full output closure inherited, consumed sources freshly hashed. All biological aims remain open.')
    atomic(a.receipt, result)
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__':main()
