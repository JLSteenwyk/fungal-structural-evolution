#!/usr/bin/env python3
"""Scalar chain screening only; never sufficient for ancestral qualification."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import arviz as az
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def diagnose(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[0] != 4:
        raise ValueError('Exactly four chains required')
    if values.shape[1] < 20:
        return dict(status='insufficient_retained_draws', draws_per_chain=values.shape[1])
    if not np.isfinite(values).all():
        return dict(status='nonfinite_draws', draws_per_chain=values.shape[1])
    if np.any(np.ptp(values, axis=1) == 0):
        return dict(status='constant_chain_requires_review', draws_per_chain=values.shape[1])
    raw = dict(rhat=float(az.rhat(values, method='rank')),
               bulk_ess=float(az.ess(values, method='bulk')),
               tail_ess=float(az.ess(values, method='tail')),
               mean_mcse=float(az.mcse(values, method='mean')))
    finite = all(np.isfinite(value) for value in raw.values())
    passed = finite and raw['rhat'] < 1.01 and min(raw['bulk_ess'], raw['tail_ess']) >= 400
    return dict(status='passes_scalar_screen_only' if passed else 'scalar_mixing_requires_review',
                draws_per_chain=values.shape[1],
                **{k: v if np.isfinite(v) else None for k, v in raw.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    chains = manifest['chains']
    assert len(chains) == len({r['chain_id'] for r in chains}) == len({r['seed'] for r in chains}) == 4
    assert len({r['model_input_identity'] for r in chains}) == 1
    assert len({str(Path(r['log']).resolve()) for r in chains}) == 4
    expected = manifest['expected_iterations']
    assert expected == sorted(set(expected)) and expected[0] == 0
    burnin = manifest['discard_through_iteration']
    assert 0 <= burnin < expected[-1]
    variables = manifest['variables']
    assert variables and len(set(variables)) == len(variables)
    values = {v: [] for v in variables}
    pins = {str(args.manifest): sha(args.manifest), str(Path(__file__)): sha(__file__)}
    for chain in chains:
        assert sha(chain['log']) == chain['log_sha256']
        pins[chain['log']] = chain['log_sha256']
        with open(chain['log']) as handle:
            records = list(csv.DictReader(handle, delimiter='\t'))
        assert [int(r['iter']) for r in records] == expected
        retained = [r for r in records if int(r['iter']) > burnin]
        for variable in variables:
            values[variable].append([float(r[variable]) for r in retained])
    results = {variable: diagnose(draws) for variable, draws in values.items()}
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    payload = dict(status='scalar_diagnostics_complete_not_posterior_qualification',
                   arviz_version=az.__version__, numpy_version=np.__version__,
                   discard_through_iteration=burnin, variables=results, pins=pins,
                   thresholds=dict(rhat_strict_upper=1.01, bulk_ess_minimum=400, tail_ess_minimum=400),
                   scope='Declared scalar variables only. Input identity and seeds are manifest claims '
                   'requiring independent provenance audit. No alignment/ancestral-state mixing, '
                   'root uncertainty, burn-in sensitivity or full posterior qualification.')
    (out / 'diagnostics.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
    print(json.dumps({v: r['status'] for v, r in results.items()}))


if __name__ == '__main__':
    main()
