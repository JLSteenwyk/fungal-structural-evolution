#!/usr/bin/env python3
"""Exercise V9 pre-dispatch metadata on three actual small SciPy SVD calls."""
import argparse
from datetime import datetime, timezone
import inspect
import json
from pathlib import Path

import numpy as np
from scipy import linalg
import scipy.linalg._flapack as flapack

from ancestral_chain_attempt import sha, write_json
import diagnose_retained_timing_watched_traced_factors_v9 as observer
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args(); assert not a.receipt.exists(); a.output.mkdir(exist_ok=False)
    pins = {}; modules = project_sources(pins, [Path(__file__)])
    bind(pins, inspect.getsourcefile(inspect.unwrap(linalg.svd)))
    bind(pins, flapack.__file__)
    calls = []
    real_probe = observer.probe
    try:
        for case in ['c_order', 'f_order', 'source_alias']:
            operand = np.array(np.arange(21, dtype=float).reshape(7, 3)+0.25,
                order='F' if case == 'f_order' else 'C')
            source = dict(factors={'artificial_a': operand if case == 'source_alias' else operand.copy(order='C'),
                'artificial_b': np.ones((7, 2), dtype=float)})
            before = {name: value.tobytes() for name, value in source['factors'].items()}
            def fixture_probe(source, matrix):
                # Only entry dispatch is artificial; this is an actual SciPy/LAPACK call.
                left, values, right = linalg.svd(matrix, full_matrices=False,
                    lapack_driver='gesvd', check_finite=True, overwrite_a=False)
                assert np.allclose((left*values)@right, matrix, rtol=1e-12, atol=1e-12)
                return dict(status='actual_small_svd_fixture', scientific_eligibility=False)
            observer.probe = fixture_probe
            events_path = a.output/(case+'.jsonl')
            with events_path.open('x') as handle:
                outcome = observer.traced_probe(source, (operand,), handle, 1, 'artificial-'+case)
            events = [json.loads(line) for line in events_path.read_text().splitlines()]
            dispatches = [e['native_svd_dispatch'] for e in events if 'native_svd_dispatch' in e]
            assert dispatches; dispatch = dispatches[0]
            assert all(item == dispatch for item in dispatches)
            assert dispatch['operand_shape'] == [7, 3] and dispatch['operand_dtype'] == 'float64'
            assert dispatch['lapack_driver'] == 'gesvd' and 'dgesvd' in dispatch['lapack_callable']
            assert dispatch['workspace_elements'] > 0 and dispatch['compute_uv']
            assert not dispatch['full_matrices'] and not dispatch['overwrite_a']
            assert dispatch['c_contiguous'] == (case != 'f_order')
            assert dispatch['f_contiguous'] == (case == 'f_order')
            assert dispatch['operand_pointer'] == operand.__array_interface__['data'][0]
            assert dispatch['operand_shares_source_factor'] == {
                'artificial_a': case == 'source_alias', 'artificial_b': False}
            assert before == {name: value.tobytes() for name, value in source['factors'].items()}
            assert all(not e['changed'] for e in events)
            calls.append(dict(case=case, actual_native_svd_calls=1, dispatch_line_events=len(dispatches), dispatch=dispatch,
                inputs_preserved=True, line_observer_events=outcome['line_observer_events']))
    finally:
        observer.probe = real_probe
    assert observer.probe is real_probe
    for path in a.output.iterdir(): bind(pins, path)
    verify(pins)
    result = dict(status='passed_actual_native_v9_dispatch_metadata_controls_v2',
        checked_utc=datetime.now(timezone.utc).isoformat(), actual_native_svd_calls=3,
        artificial_probe_entry=True, original_probe_and_guard_runs=0,
        cases=calls, transitive_project_source_modules=len(modules), source_hashes=pins,
        original_jobs_restarted=False, production_tolerance_changed=False,
        biological_fits=0, scientific_eligibility=False, gpu=False,
        scope='Three actual small SciPy gesvd calls: C order, F order and direct source-factor alias. '
              'Repeated CPython line events retain identical pre-dispatch metadata; they are not counted as extra native calls. Shape/type/workspace/flags/pointer/alias metadata checked. Probe entry alone '
              'is substituted in this private software process and restored; original math/guards are '
              'not exercised here. No original factors, mutation attribution, repair or biology.')
    write_json(a.receipt, result)
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__': main()
