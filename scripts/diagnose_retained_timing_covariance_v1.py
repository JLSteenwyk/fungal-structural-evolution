#!/usr/bin/env python3
"""Preserve exact failing-cohort guards and independent covariance comparisons.

This diagnostic records guard rejection as evidence, never as fit acceptance.
No inherited audit, production tolerance or original output is modified.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import time

import numpy as np
from scipy import sparse

from ancestral_chain_attempt import sha
from covariance_basis_context import ComponentKernelProducts
from covariance_basis_independent import IndependentKernelProducts
from full_covariance_qualification_sources import folded_operators
from full_expanded_model_design_sources import array_digest, digest
from full_retained_shared_entity_fit_sources import cohorts, operators_for
from full_retained_shared_entity_timing import groups
from prepare_full_retained_shared_entity_timing import sources
from reference_measurement_union_sources import bind, verify
from retained_shared_entity_candidate import backend_guard, validate_source
from shared_entity_likelihood import SharedEntityLikelihood


def save(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def comparison(a, b, names, error_a=None, error_b=None):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    delta = abs(a - b); tolerance = 2e-8 + 3e-9 * abs(b)
    mask = delta > tolerance
    entries = [dict(row=int(i), column=int(j), row_name=names[i], column_name=names[j],
                    actual=float(a[i,j]), reference=float(b[i,j]),
                    absolute_difference=float(delta[i,j]), tolerance=float(tolerance[i,j]))
               for i, j in np.argwhere(mask)]
    result = dict(rtol=3e-9, atol=2e-8, mismatched_entries=entries,
                  passes_original_allclose=not bool(mask.any()),
                  maximum_absolute_difference=float(delta.max()),
                  maximum_relative_difference=float(np.max(delta / np.maximum(abs(b), np.finfo(float).tiny))))
    if error_a is not None:
        envelope = np.asarray(error_a) + np.asarray(error_b)
        rounding = 8 * np.finfo(float).eps * np.maximum(abs(a), abs(b))
        result['preserved_envelope_excess_entries'] = np.argwhere(delta > envelope + rounding).tolist()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    plan = json.loads(args.plan.read_text()); verify(plan['pins'])
    root = Path(plan['output']); root.mkdir(exist_ok=False)
    save(root / 'stage_plan.json', plan)
    timing_path = Path(plan['original_timing_plan'])
    timing = json.loads(timing_path.read_text())
    print('diagnostic_source_load_started', flush=True)
    fit, source, bindings = sources(timing, timing_path)
    for path, checksum in plan['pins'].items(): bind(bindings, path, checksum)
    bind(bindings, args.plan)
    print('diagnostic_source_load_complete', len(bindings), flush=True)
    target = plan['cohort_id']
    index = next(i for i, c in enumerate(source['cohorts'], 1) if c['cohort_id'] == target)
    assert index == 11
    for number, (cohort, rows, entries) in enumerate(cohorts(source, fit), 1):
        print('diagnostic_verified_cohort_prefix', number, flush=True)
        if cohort['cohort_id'] == target: break
    else: raise AssertionError('Original failing cohort absent')
    assert len(rows) == plan['records'] == 22881
    census, selected, counts = groups(source, fit, cohort, entries)
    assert len(census) == 1200 and len(selected) == 40
    save(root / 'candidate_census.json', census)
    save(root / 'cohort.json', cohort)
    np.save(root / 'case_rows.npy', rows)
    np.save(root / 'component_labels.npy', source['labels'][rows])
    save(root / 'case_ids.json', [source['ids'][i] for i in rows])
    operators = {}; full_banks = {}; factors = {}; raw_contexts = {}; latent_contexts = {}
    for mode in fit['loading_modes']:
        rep = next(r for r in selected.values() if r['key'][0] == mode)
        operators[mode] = operators_for(source, rows, rep['source_audit'])
        for name, z in operators[mode].items(): sparse.save_npz(root / (mode + '-' + name + '.npz'), z)
        full_banks[mode] = ComponentKernelProducts(source['labels'][rows], folded_operators(source, rows, mode), np.ones(len(rows)))
        retained_bank = ComponentKernelProducts(source['labels'][rows], operators[mode], np.ones(len(rows)))
        latent_bank = IndependentKernelProducts(source['labels'][rows], operators[mode])
        for tree in fit['trees']:
            if tree not in factors:
                factors[tree] = source['factors'][tree][rows]
                np.save(root / (tree + '-factor.npy'), factors[tree])
            raw_contexts[mode, tree] = retained_bank.tree(factors[tree])
            latent_contexts[mode, tree] = latent_bank.tree(factors[tree])
    folder = root / 'groups'; folder.mkdir()
    results = []; first_failure = None
    for group_id, rep in sorted(selected.items()):
        identity = rep['identity']; mode, tree, method, outcome = rep['key']
        x = rep['matrix']; y = rep['response']; audit = rep['source_audit']
        validate_source(identity, audit, rep['original_audit'], rep['certificate'])
        xhash = array_digest(x, '<f8'); fhash = array_digest(factors[tree], '<f8')
        primary = SharedEntityLikelihood(source['labels'][rows], operators[mode], factors[tree], np.ones(len(rows)), x, y)
        assert array_digest(x, '<f8') == xhash and array_digest(factors[tree], '<f8') == fhash
        del primary
        # Preserve the exact original guard call and its exception. Catching it
        # here makes a diagnostic record, never a passing production candidate.
        try:
            backend_guard(source, rows, x, operators[mode], audit)
            guard = dict(status='passed_unchanged_original_guard')
        except (AssertionError, ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            guard = dict(status='rejected_by_unchanged_original_guard', error_type=type(error).__name__, error_message=str(error))
            if first_failure is None: first_failure = group_id
        fresh = raw_contexts[mode, tree].audit(x)
        latent = latent_contexts[mode, tree].project(x)
        full = full_banks[mode].tree(factors[tree]).audit(x)
        names = audit['retained_kernel_names']; subset = [full['kernel_names'].index(n) for n in names]
        inherited = audit['numerical_audit']
        comparisons = {}
        for key, error, latent_key, latent_error in [
            ('raw_gram', 'raw_roundoff_envelope', 'raw', 'raw_error'),
            ('projected_gram', 'projected_roundoff_envelope', 'projected', 'projected_error')]:
            principal = np.asarray(full[key])[np.ix_(subset, subset)]
            principal_error = np.asarray(full[error])[np.ix_(subset, subset)]
            comparisons[key] = dict(
                fresh_vs_inherited=comparison(fresh[key], inherited[key], names, fresh[error], inherited[error]),
                latent_vs_fresh=comparison(latent[latent_key], fresh[key], names, latent[latent_error], fresh[error]),
                latent_vs_inherited=comparison(latent[latent_key], inherited[key], names, latent[latent_error], inherited[error]),
                full_principal_vs_inherited=comparison(principal, inherited[key], names, principal_error, inherited[error]),
                full_principal_vs_fresh=comparison(principal, fresh[key], names, principal_error, fresh[error]))
        record = dict(group_id=group_id, representative=identity, eligible_candidates=counts[group_id],
                      selection_rank=list(rep['rank']), guard=guard, comparisons=comparisons,
                      design_sha256=xhash, factor_sha256=fhash, factor_dtype=str(factors[tree].dtype),
                      primary_constructor_preserves_inputs=True, scientific_eligibility=False)
        detail = dict(**record, fresh=fresh, inherited=inherited, latent={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in latent.items()},
                      full_original_seven_kernel_recomputation=full,
                      original_source_audit=rep['original_audit'], retained_source_audit=audit, certificate=rep['certificate'])
        np.savez_compressed(folder / (group_id + '.npz'), design=x, response=y)
        save(folder / (group_id + '.json'), detail)
        results.append(record)
        print('diagnostic_group', len(results), '/', len(selected), group_id, guard['status'], flush=True)
    verify(bindings)
    archive = root / 'source_hashes.json.gz'
    with gzip.open(archive, 'xt') as handle: json.dump(bindings, handle, sort_keys=True)
    artifacts = {str(p):sha(p) for p in root.rglob('*') if p.is_file()}
    receipt = dict(status='completed_exact_failed_cohort_covariance_diagnostic_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(), plan=str(args.plan), plan_sha256=sha(args.plan),
        cohort_id=target, records=len(rows), verified_cohort_prefix=index,
        candidate_census_rows=len(census), selected_groups=len(results), first_sorted_guard_failure=first_failure,
        guard_status_counts=dict(Counter(r['guard']['status'] for r in results)), groups=results,
        full_source_bindings_freshly_verified=len(bindings), full_source_hash_archive=str(archive),
        full_source_hash_archive_sha256=sha(archive), artifacts=artifacts,
        elapsed_seconds=time.monotonic()-started, scientific_eligibility=False,
        production_tolerance_changed=False, original_jobs_restarted=False, biological_fits=0, gpu=False,
        scope='Exact original cohort11, original selection of40representatives from1200candidates, '
              'source validation and primary construction followed by unchanged original guard. '
              'All rejections retained. Separate retained/component, independent latent and full '
              'seven-kernel principal-submatrix arithmetic comparisons. Diagnostic only: no '
              'promotion of inherited source, biological fit, tolerance relaxation or original run restart.')
    save(args.receipt, receipt)
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['groups','artifacts']}, indent=2), flush=True)


if __name__ == '__main__': main()
