#!/usr/bin/env python3
"""Cost and bind the first full independent-chain sampling horizon."""
import csv
import json
import math
from pathlib import Path
import shutil
from ancestral_chain_attempt import sha, write_json


def main():
    inputs = Path('results/ancestral/baliphy-independent-chain-inputs-20260927-v1')
    resources = Path('results/ancestral/baliphy-resources-20260927-final-v1')
    pins = {}
    for root in [inputs, resources]:
        receipt = json.loads((root / 'receipt.json').read_text())
        pins[str(root / 'receipt.json')] = sha(root / 'receipt.json')
        for path, digest in receipt['pins'].items():
            assert sha(path) == digest, path
        for path, digest in receipt['artifacts'].items():
            assert sha(root / path) == digest
            pins[str(root / path)] = digest
    summary = json.loads((resources / 'receipt.json').read_text())
    assert summary['audit_dispositions'] == {'all_saved_samples_and_candidate_nodes_verified': 324}
    groups = {r['group']: r for r in csv.DictReader((resources / 'effective_group_resources.tsv').open(), delimiter='\t')}
    chains = json.loads((inputs / 'chain_inputs.json').read_text())
    assert len(chains) == 1620 and len(groups) == 135
    binary = str(Path(json.loads(Path('metadata/baliphy_prior_initialization_plan_20260927.json').read_text())['binary']).resolve())
    prlimit = str(Path(shutil.which('prlimit')).resolve())
    iterations, workers = 1000, 16
    jobs = []
    for chain in sorted(chains, key=lambda r: (r['chain'], r['prior_label'], -float(groups[r['effective_input_group']]['median_measured_seconds']), r['chain_id'])):
        seconds = float(groups[chain['effective_input_group']]['median_measured_seconds'])
        # Four-fold timing allowance is a planning bound, not a calibrated interval.
        timeout = max(3600, math.ceil(seconds * iterations / 20 * 4))
        paths = [binary, prlimit, chain['program'], chain['tree'], chain['alignment']]
        config = dict(command=[prlimit, '--as=' + str(12 * 2**30), '--fsize=' + str(2 * 2**30), '--',
            binary, '--seed', str(chain['seed']), 'run', chain['program'],
            '--iterations', str(iterations), '--name', 'independent-chain'],
            seed=chain['seed'], model_input_identity=chain['effective_input_group'] + '-' + chain['prior_label'],
            timeout_seconds=timeout, pins={str(Path(p).resolve()): sha(p) for p in paths})
        jobs.append(dict(chain=chain, config=config))
    jobs_path = inputs / 'initial_horizon_jobs.json'
    assert not jobs_path.exists()
    write_json(jobs_path, jobs)
    pins[str(jobs_path)] = sha(jobs_path)
    mapping = 'results/ancestral/case-local-trees-20260927-v1/ancestral_node_mapping.tsv'
    paths = [mapping, __file__, 'scripts/run_independent_baliphy_chains.py',
        'scripts/readback_independent_baliphy_chain.py', 'scripts/ancestral_chain_attempt.py',
        'scripts/audit_baliphy_sample_mapping.py', 'scripts/prepare_case_ancestral_neighborhoods.py',
        'metadata/baliphy_sample_mapping_final_readback_completed_20260927.json',
        'metadata/baliphy_resource_final_completed_20260927.json',
        'metadata/baliphy_attempt_recovery_completed_20260927.json',
        'metadata/baliphy_independent_readback_fixture_check_20260927.json',
        'scripts/ancestral_chain_diagnostics.py', 'environments/ancestral-diagnostics-20260927.lock.txt']
    pins.update({str(p): sha(p) for p in paths})
    baseline_hours = sum(float(g['median_measured_seconds']) for g in groups.values()) * iterations / 20 * 12 / 3600
    plan = dict(unit='fungal-baliphy-independent-chains-20260927.service', jobs=str(jobs_path),
        output='results/ancestral/baliphy-independent-chains-20260927-v1',
        workers=workers, iterations=iterations, memory_bytes=192 * 2**30,
        minimum_free_disk_bytes=3072 * 2**30, mapping=mapping, pins=pins,
        completed_units=['fungal-baliphy-sample-mapping-20260927.service',
                         'fungal-baliphy-sample-mapping-final-audit-20260927.service',
                         'fungal-baliphy-resource-summary-20260927.service'],
        resources=dict(cpus=workers, memory_gib=192, swap_gib=0, per_process_address_space_gib=12,
            per_file_limit_gib=2, output_planning_gib=256, paid_cost=0, gpu=False,
            fixed_parameter_linear_worker_hours=baseline_hours,
            scheduling_days_at_linear_baseline=baseline_hours / workers / 24,
            sensitivity_factor_range=[.5, 4],
            timeout_sum_worker_hours=sum(j['config']['timeout_seconds'] for j in jobs) / 3600,
            caveat='Startup-inclusive 20-iteration fixed-parameter runs scaled to free-parameter chains. '
                   'Not a reliable ETA or bound on convergence; timeouts remain unresolved outcomes.'),
        diagnostics=dict(burn_in_fractions=[.25, .5], rhat_strict_upper=1.01,
            bulk_ess_minimum=400, tail_ess_minimum=400,
            policy='Evaluate both cutoffs and all stochastic scalars; separately assess saved alignment, '
                   'candidate length and ancestral-state mixing. Failed screens require longer fresh-chain '
                   'horizons or sampler review, retaining prior attempts. Do not concatenate or declare '
                   'convergence from this horizon. Root and prior sensitivity remain explicit.'),
        scope='All135 effective inputs x3 priors x4 seeds; all324 aliases retained. '
              'First sampling horizon, not a pilot, converged posterior, or ancestral structural ensemble.')
    path = Path('metadata/baliphy_independent_chain_plan_20260927.json')
    assert not path.exists()
    write_json(path, plan)
    print(json.dumps(plan['resources'], indent=2))


if __name__ == '__main__':
    main()
