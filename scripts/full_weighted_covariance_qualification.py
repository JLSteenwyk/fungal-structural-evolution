#!/usr/bin/env python3
"""Fresh complete-grid weighted raw/REML qualification and latent readback.

Every original setting is retained. Segments bound individual output files;
they do not select cohorts, controls, trees, settings or numerical cases.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import itertools
import json
import os
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from full_exact_covariance_sources import closed_subset
from full_expanded_model_design_sources import AXES, DEGREES, SETTING_FIELDS, array_digest, digest
from full_expanded_model_input_sources import ORDERS
from full_weighted_covariance_sources_v2 import MODES, POLICIES, basis, cohort_controls, design_matrix, jsonl, load, retained_operators
from independent_positive_diagonal_basis_context import IndependentPositiveDiagonalBasisContext
from positive_diagonal_basis_context import PositiveDiagonalBasisContext
from readback_full_covariance_qualification import numeric
from reference_measurement_union_sources import bind, verify

SCHEMA = 'full-four-control-covariance-numerical-qualification-v1'
PRODUCER = 'complete_full_four_control_covariance_numerical_qualification_pending_readback_v1'
READER = 'passed_full_four_control_covariance_numerical_latent_readback_v1'
LINK_FIELDS = ['source_setting_ordinal', 'source_setting_sha256', *SETTING_FIELDS,
    'loading_mode', 'tree', 'control_policy', 'audit_id', 'covariance_disposition', 'combined_disposition']
SUMMARY = ['logical_cases', 'cohorts', 'designs', 'settings', 'numerical_audit_rows', 'setting_audit_links',
    'audit_status_counts', 'link_status_counts', 'policy_audit_status_counts', 'policies', 'trees', 'loading_modes']


def write(path, value):
    with path.open('x') as f: json.dump(value, f, sort_keys=True, allow_nan=False); f.write('\n')


def identity(contract, did, mode, tree, policy):
    return digest([SCHEMA, contract, did, mode, tree, policy])


def disposition(value):
    return 'numerically_qualified_exact_four_control_covariance_basis' if all(
        value[k]['disposition'] == 'numerically_independent_covariance_bases'
        for k in ['raw_diagnostics', 'reml_diagnostics']) else 'exact_four_control_covariance_basis_requires_review'


def source_inputs(plan, path):
    source, bindings = load(plan, path)
    cp = Path(plan['source_census_plan']); spec = json.loads(cp.read_text()); root = Path(spec['output'])
    completion = closed_subset(plan['source_census_completion'],
        'complete_verified_full_four_control_covariance_source_census_v1',
        [cp, root / 'receipt.json', root / 'readback.json', root / 'cohort_source_census.jsonl.gz'], bindings)
    assert completion['cohorts'] == plan['expected']['cohorts']
    for field, original in [('designs', 'designs'), ('settings', 'settings'),
        ('case_row_occurrences', 'case_row_occurrences')]: assert completion[field] == plan['expected'][original]
    assert completion['expected_numerical_audit_rows'] == plan['expected']['numerical_audit_rows']
    assert completion['expected_setting_audit_links'] == plan['expected']['setting_audit_links']
    for k in ['operator', 'design', 'inputs', 'covariance', 'exact', 'cone', 'weights']:
        assert spec[k + '_plan'] == plan[k + '_plan'] and spec[k + '_completion'] == plan[k + '_completion']
    assert spec['trees'] == plan['trees']
    source['contract'] = digest(dict(schema=SCHEMA, source_contract=source['contract'],
        source_census_completion_sha256=sha(plan['source_census_completion']),
        source_hashes=bindings, numerical_audits_inherited=False))
    verify(bindings)
    return source, bindings


def runtime(plan):
    group = next(l[3:] for l in Path('/proc/self/cgroup').read_text().splitlines() if l.startswith('0::'))
    root = Path('/sys/fs/cgroup') / group.lstrip('/'); r = plan['resources']
    assert {k:(root / k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']} == {
        'cpu.max': str(r['cpus'] * 100000) + ' 100000', 'memory.max': str(r['memory_gib'] * 2**30), 'memory.swap.max': '0'}
    assert all(os.environ.get(k) == '1' for k in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'])


def settings_by_cohort(source):
    groups = {c['cohort_id']: [] for c in source['cohorts']}; settings = 0
    with gzip.open(source['root'] / 'model_settings.tsv.gz', 'rt') as f:
        reader = csv.reader(f, delimiter='\t'); assert next(reader) == SETTING_FIELDS
        index = SETTING_FIELDS.index('cohort_id')
        for ordinal, row in enumerate(reader):
            assert len(row) == len(SETTING_FIELDS)
            groups[row[index]].append((ordinal, row)); settings += 1
    return groups, settings


def expected_record(source, cohort, design, route, diagonal, control, mode, tree, policy):
    return dict(schema=SCHEMA, audit_id=identity(source['contract'], design['design_id'], mode, tree, policy),
        source_contract=source['contract'], cohort_id=cohort['cohort_id'], design_id=design['design_id'],
        loading_mode=mode, tree=tree, control_policy=policy, records=cohort['records'],
        source_design_disposition=design['disposition'], raw_design_sha256=design['raw_design_sha256'],
        active_column_indices=design['active_column_indices'], cohort_rows_sha256=cohort['case_rows_sha256'],
        ordered_case_ids_sha256=cohort['ordered_case_ids_sha256'], diagonal_sha256=array_digest(diagonal, '<f8'),
        control_record_sha256=control['record_sha256'], retained_kernel_names=route['names'],
        exact_certificate_sha256=route['certificate_sha256'], basis_route=route['route'],
        residual_diagonal_is_exact_uniform_one=route['exact_uniform_one'],
        scientific_eligibility=False, nonuniform_weighting_accepted=False,
        component_variance_attribution_accepted=False, covariance_model_selected=False)


def failure_capture(root, expected, saved, reference, selected, diagonal, design, factor, error):
    directory = root / 'failures' / expected['audit_id']; directory.mkdir(parents=True, exist_ok=False)
    write(directory / 'failure.json', dict(expected_identity=expected, saved=saved,
        reference={k:v.tolist() if isinstance(v, np.ndarray) else v for k,v in reference.items()} if reference else None,
        error_type=type(error).__name__, error_message=str(error), scientific_eligibility=False))
    with (directory / 'original_numeric_inputs.npz').open('xb') as f:
        np.savez_compressed(f, case_rows=selected, diagonal=diagonal, design=design, species_factor=factor)


def run(path, reader=False):
    plan = json.loads(path.read_text()); runtime(plan); source, bindings = source_inputs(plan, path)
    original = dict(bindings); root = Path(plan['output']); receipt_name = 'readback.json' if reader else 'receipt.json'
    assert not (root / receipt_name).exists()
    if reader:
        receipt = json.loads((root / 'receipt.json').read_text())
        assert receipt['status'] == PRODUCER and receipt['source_contract'] == source['contract']
        assert receipt['plan_sha256'] == sha(path) and receipt['source_hashes'] == original
        assert receipt['scientific_eligibility'] is False
        bind(bindings, root / 'receipt.json')
        for n,h in receipt['artifacts'].items(): bind(bindings, root / n, h)
        manifest = json.loads((root / 'cohort_manifest.json').read_text())
        assert [r['cohort_id'] for r in manifest] == [c['cohort_id'] for c in source['cohorts']]
        assert json.loads((root / 'stage_plan.json').read_text()) == dict(schema=SCHEMA, plan_sha256=sha(path), source_contract=source['contract'])
    else:
        root.mkdir(parents=True, exist_ok=False); (root / 'cohorts').mkdir()
        write(root / 'stage_plan.json', dict(schema=SCHEMA, plan_sha256=sha(path), source_contract=source['contract']))
        manifest = []
    setting_groups, settings_n = settings_by_cohort(source)
    assert settings_n == plan['expected']['settings']
    designs = jsonl(source['root'] / 'unique_designs.jsonl')
    counts = Counter(); link_counts = Counter(); policy_counts = Counter(); artifacts = {}
    audit_n = link_n = design_n = conservative_reviews = 0; maximum_error = 0.
    for ci, (cohort, control) in enumerate(zip(source['cohorts'], source['controls'])):
        selected, diagonals, control_record = cohort_controls(source, cohort, control)
        cohort_designs = [next(designs) for _ in range(30)]
        assert {(d['order_contrast'], d['sequence_axis'], d['degree']) for d in cohort_designs} == set(itertools.product(ORDERS, AXES, DEGREES))
        matrices = {d['design_id']: design_matrix(source, cohort, selected, d) for d in cohort_designs}
        design_n += len(cohort_designs); statuses = {}; cohort_counts = Counter(); cohort_link_counts = Counter()
        ap = root / 'cohorts' / (cohort['cohort_id'] + '.audits.jsonl.gz')
        lp = root / 'cohorts' / (cohort['cohort_id'] + '.links.tsv.gz')
        exports = jsonl(ap) if reader else None
        out = None if reader else gzip.open(ap, 'wt')
        try:
            for mode, tree in itertools.product(MODES, plan['trees']):
                routes = [basis(source, cohort, mode, d) for d in diagonals]
                contexts = {}; operators_by_basis = {}
                factor = source['factors'][tree][selected]
                for route in routes:
                    key = tuple(route['names'])
                    if key not in contexts:
                        ops = retained_operators(source, selected, mode, route['names']); operators_by_basis[key] = ops
                        cls = IndependentPositiveDiagonalBasisContext if reader else PositiveDiagonalBasisContext
                        contexts[key] = cls(source['labels'][selected], ops, factor)
                for d in cohort_designs:
                    x = matrices[d['design_id']]; projected = {}; context_errors = {}
                    if d['disposition'] == 'full_rank_design':
                        for key, context in contexts.items():
                            try: projected[key] = context.design(x)
                            except (ValueError, ArithmeticError) as e: context_errors[key] = e
                    for j, policy in enumerate(POLICIES):
                        route = routes[j]; key = tuple(route['names']); diagonal = diagonals[j]
                        expected = expected_record(source, cohort, d, route, diagonal, control, mode, tree, policy)
                        if reader:
                            saved = next(exports); reference = None
                            try:
                                assert all(saved[k] == v for k,v in expected.items())
                                if d['disposition'] != 'full_rank_design':
                                    assert saved['disposition'] == d['disposition'] and saved['numerical_audit'] is None
                                elif saved['numerical_audit'] is None:
                                    assert saved['disposition'] == 'numerical_covariance_qualification_requires_review'
                                    primary = PositiveDiagonalBasisContext(source['labels'][selected], operators_by_basis[key], factor)
                                    try: primary.design(x).audit(diagonal)
                                    except (ValueError, ArithmeticError) as e:
                                        assert saved['error_type'] == type(e).__name__ and saved['error_message'] == str(e)
                                    else: raise AssertionError('False primary numerical failure')
                                else:
                                    assert key in projected, 'Independent design requires review'
                                    reference = projected[key].audit(diagonal); value = saved['numerical_audit']
                                    assert value['fixed_effect_columns'] == x.shape[1] and value['residual_dimension'] == len(selected) - x.shape[1]
                                    assert value['family_components'] == len(np.unique(source['labels'][selected]))
                                    assert value['exactly_zero_incidence_names'] == [k for k,z in operators_by_basis[key].items() if not z.nnz]
                                    # Frozen verifier uses supplied kernels and fresh D bounds;
                                    # its historical uniform disposition name is never exported.
                                    review, error = numeric(value, reference, route['names'], len(selected))
                                    conservative_reviews += int(review); maximum_error = max(maximum_error, error)
                                    assert saved['disposition'] == disposition(value)
                            except (AssertionError, ValueError, ArithmeticError, KeyError) as e:
                                failure_capture(root, expected, saved, reference, selected, diagonal, x, factor, e)
                                raise
                        else:
                            saved = dict(expected, numerical_audit=None)
                            if d['disposition'] != 'full_rank_design': saved['disposition'] = d['disposition']
                            else:
                                try:
                                    if key in context_errors: raise context_errors[key]
                                    saved['numerical_audit'] = projected[key].audit(diagonal)
                                    saved['disposition'] = disposition(saved['numerical_audit'])
                                except (ValueError, ArithmeticError) as e:
                                    saved.update(disposition='numerical_covariance_qualification_requires_review',
                                        error_type=type(e).__name__, error_message=str(e))
                            out.write(json.dumps(saved, sort_keys=True, allow_nan=False) + '\n')
                        status = saved['disposition']; aid = expected['audit_id']
                        assert (d['design_id'], mode, tree, policy) not in statuses
                        statuses[d['design_id'], mode, tree, policy] = aid, status
                        counts[status] += 1; cohort_counts[status] += 1; policy_counts[policy + ':' + status] += 1; audit_n += 1
        finally:
            if out is not None: out.close()
        if reader: assert next(exports, None) is None
        assert len(statuses) == 1200
        def expected_links():
            for ordinal, row in setting_groups[cohort['cohort_id']]:
                setting = dict(zip(SETTING_FIELDS, row)); assert setting['cohort_id'] == cohort['cohort_id']
                for mode, tree, policy in itertools.product(MODES, plan['trees'], POLICIES):
                    aid, status = statuses[setting['design_id'], mode, tree, policy]
                    combined = setting['disposition'] if setting['disposition'] != 'ready_for_working_covariance_fit' else status
                    yield [str(ordinal), digest(row), *row, mode, tree, policy, aid, status, combined]
        with gzip.open(lp, 'rt' if reader else 'wt') as f:
            if reader:
                links = csv.reader(f, delimiter='\t'); assert next(links) == LINK_FIELDS
            else:
                links = csv.writer(f, delimiter='\t', lineterminator='\n'); links.writerow(LINK_FIELDS)
            for row in expected_links():
                if reader: assert next(links) == row
                else: links.writerow(row)
                link_counts[row[-1]] += 1; cohort_link_counts[row[-1]] += 1; link_n += 1
            if reader: assert next(links, None) is None
        entry = dict(cohort_id=cohort['cohort_id'], audits=1200, links=len(setting_groups[cohort['cohort_id']]) * 40,
            audit_file=str(ap.relative_to(root)), audit_sha256=sha(ap), link_file=str(lp.relative_to(root)), link_sha256=sha(lp),
            audit_status_counts=dict(cohort_counts), link_status_counts=dict(cohort_link_counts))
        if reader: assert manifest[ci] == entry
        else: manifest.append(entry)
        artifacts[entry['audit_file']] = entry['audit_sha256']; artifacts[entry['link_file']] = entry['link_sha256']
        if (ci + 1) % 10 == 0: print('independently_checked_weighted_covariance_cohorts' if reader else 'fresh_weighted_covariance_cohorts', ci + 1, flush=True)
    assert next(designs, None) is None
    mp = root / 'cohort_manifest.json'
    if not reader: write(mp, manifest)
    artifacts['cohort_manifest.json'] = sha(mp); artifacts['stage_plan.json'] = sha(root / 'stage_plan.json')
    summary = dict(logical_cases=len(source['ids']), cohorts=len(source['cohorts']), designs=design_n, settings=settings_n,
        numerical_audit_rows=audit_n, setting_audit_links=link_n, audit_status_counts=dict(counts), link_status_counts=dict(link_counts),
        policy_audit_status_counts=dict(policy_counts), policies=POLICIES, trees=plan['trees'], loading_modes=MODES)
    assert design_n == plan['expected']['designs'] and audit_n == plan['expected']['numerical_audit_rows']
    assert link_n == plan['expected']['setting_audit_links']
    if reader:
        assert all(receipt[k] == v for k,v in summary.items()) and receipt['artifacts'] == artifacts
    verify(bindings)
    result = dict(status=READER if reader else PRODUCER, checked_utc=datetime.now(timezone.utc).isoformat(),
        plan_sha256=sha(path), source_contract=source['contract'], **summary, source_hashes=bindings, artifacts=artifacts,
        working_model_fits_computed=0, scientific_eligibility=False, nonuniform_weighting_accepted=False,
        scope=plan['scope'])
    if reader: result.update(producer_receipt_sha256=sha(root / 'receipt.json'),
        conservative_diagnostic_disagreements=conservative_reviews, maximum_projected_gram_difference=maximum_error)
    write(root / receipt_name, result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes', 'artifacts', 'scope']}), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--reader', action='store_true'); a = p.parse_args(); run(a.plan, a.reader)
