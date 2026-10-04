"""Actual native guard checks and private full-export mutation negatives."""
from contextlib import nullcontext
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy import sparse

from check_full_shared_entity_fits_v2 import settings
from check_retained_shared_entity_candidates import inputs
from check_full_weighted_shared_entity_fits_v2 import synthetic_candidate, synthetic_reader
from full_weighted_covariance_qualification_parallel_v1 import arrays
from full_weighted_fit_exports import atomic
from full_weighted_shared_entity_fit_sources_parallel_v1 import cohorts as original_cohorts
from prepare_full_weighted_parallel_source_fits_v2 import run
from readback_full_weighted_parallel_source_fits_v2 import run as readback
from weighted_shared_entity_candidate import candidate
from readback_weighted_shared_entity_candidate import numeric
from weighted_fit_numeric_memory_guard_v1 import ArrayGuard, guarded_call


def actual_guarded_native(root):
    saved = json.loads(Path('data/software_audits/weighted-shared-entity-candidates-20261004-v1/actual_weighted_candidate_cases.json').read_text())
    assert len(saved) == 64
    plan = settings()
    checks = []
    for row in saved:
        source, rows, x, y, ops, _, _, _ = inputs(row['pair_exception'], row['loading_mode'])
        operators = {} if row['policy'] == 'uniform' else {'target_node':sparse.eye(len(rows), format='csr')}
        operators.update(ops)
        response = y if row['outcome'] == 'rmsd_delta' else np.tanh(y) + .15*x[:,1]
        diagonal = np.asarray(row['diagonal'])
        expected = row['candidate']
        identity = {k:v for k,v in expected.items() if k not in ['disposition','fit','numerical_attempted','error_type','error_message']}
        audit = row['audit']
        route = dict(names=audit['retained_kernel_names'], certificate_sha256=audit['exact_certificate_sha256'],
            route=audit['basis_route'], exact_uniform_one=audit['residual_diagonal_is_exact_uniform_one'])
        source_guard = ArrayGuard(source, 'actual-source')
        produced = guarded_call(candidate, source, plan, rows, identity, x, response,
            operators, audit, route, diagonal, root=root)
        assert produced == expected, 'Guarded native producer differs from frozen saved candidate'
        independent = guarded_call(numeric, source, plan, rows, identity, x, response,
            operators, audit, route, diagonal, root=root, exported=produced, reader=True)
        assert independent == row['independent']
        source_guard.check(root, identity, 'actual-complete-source-after')
        checks.append(dict(candidate_id=identity['candidate_id'], native_producer_matches_saved=True,
            native_reader_matches_saved=True, cached_source_preserved=True, scientific_eligibility=False))
    atomic(Path(root)/'actual_guarded_native_checks.json', checks)
    return len(checks)


def memory_negatives(root, parent_plan, private_copy):
    rejected = []
    for role in ['producer', 'reader']:
        for case in ['unused-factor', 'design', 'response', 'diagonal', 'rows',
                     'operator', 'design-layout', 'cohort-design', 'cached-before']:
            output = Path(root)/('memory-'+role+'-'+case)
            pp = Path(root)/('memory-'+role+'-'+case+'.plan.json')
            if role == 'reader':
                private_copy(parent_plan, pp, output)
            else:
                plan = json.loads(Path(parent_plan).read_text())
                plan['output'] = str(output)
                with pp.open('x') as handle:json.dump(plan, handle, indent=2)
            injected = []
            module = ('prepare' if role == 'producer' else 'readback')+'_full_weighted_parallel_source_fits_v2'
            original = synthetic_candidate if role == 'producer' else synthetic_reader

            def mutate(arguments):
                source, _, rows, identity, matrix, response = arguments[:6]
                offset = 0 if role == 'producer' else 1
                operators, _, _, diagonal = arguments[6+offset:10+offset]
                if case == 'unused-factor':
                    values = next(value for key,value in source['factors'].items() if key != identity['tree'])
                elif case in ['design', 'cohort-design', 'design-layout']:
                    values = matrix
                elif case == 'response':values = response
                elif case == 'diagonal':values = diagonal
                elif case == 'rows':values = rows
                elif case == 'operator':
                    values = next(value for _,value in arrays(operators, 'ops') if value.size and np.issubdtype(value.dtype, np.floating))
                else:raise AssertionError(case)
                if case == 'design-layout':
                    assert values.ndim == 2 and values.shape[0] != values.shape[1]
                    assert values.strides[0] != 0
                    values.strides = (0, *values.strides[1:])
                elif np.issubdtype(values.dtype, np.integer):values.flat[0] += 1
                else:values.flat[0] = np.nextafter(values.flat[0], np.inf)
                injected.append(case)

            def injected_call(*arguments):
                result = original(*arguments)
                if not injected:mutate(arguments)
                return result

            def cohort_stream(source, plan):
                for cohort, rows, entries in original_cohorts(source, plan):
                    if case == 'cached-before' and not injected:
                        factor = next(iter(source['factors'].values()))
                        factor.flat[0] = np.nextafter(factor.flat[0], np.inf)
                        injected.append(case)
                    yield cohort, rows, entries

            def changed_operators(source, rows, audit):
                from full_weighted_shared_entity_fit_sources_parallel_v1 import operators_for
                operators = operators_for(source, rows, audit)
                # Mutate the cached cohort matrix before its native-call
                # snapshot. Only the independent cohort boundary catches it.
                if not injected:
                    cached_matrix.flat[0] = np.nextafter(cached_matrix.flat[0], np.inf)
                    injected.append(case)
                return operators

            def stream_with_matrix(source, plan):
                nonlocal cached_matrix
                for cohort, rows, entries in cohort_stream(source, plan):
                    cached_matrix = entries[0][1]
                    yield cohort, rows, entries

            cached_matrix = None
            function = 'candidate' if role == 'producer' else 'numeric'
            with patch(module+'.'+function, side_effect=injected_call if case not in ['cohort-design','cached-before'] else original), \
                 patch(module+'.cohorts', side_effect=stream_with_matrix), \
                 patch(module+'.operators_for', side_effect=changed_operators) if case == 'cohort-design' else nullcontext():
                try:
                    if role == 'producer':run(pp)
                    else:readback(pp, output/'rejected-reader.json')
                except AssertionError:pass
                else:raise AssertionError('In-memory mutation passed '+role+' '+case)
            assert injected == [case]
            failures = list((output/'memory_failures').glob('*/failure.json'))
            assert failures, 'Mutation rejected without exact memory capture'
            for failure in failures:
                data = json.loads(failure.read_text())
                assert data['changed_arrays'] and data['scientific_eligibility'] is False
                for filename in data['captured_changed_arrays'].values():
                    with (failure.parent/filename).open('rb') as handle:np.load(handle, allow_pickle=False)
            if role == 'producer':assert not (output/'receipt.json').exists()
            else:assert not (output/'rejected-reader.json').exists()
            rejected.append(dict(role=role, case=case, failures=[str(p) for p in failures]))
    return rejected
