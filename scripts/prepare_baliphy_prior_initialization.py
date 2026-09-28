#!/usr/bin/env python3
"""Bind three explicit prior settings to all 135 effective ancestral inputs."""
import collections
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    original = Path('metadata/baliphy_initialization_plan_20260927.json')
    baseline = json.loads(original.read_text())
    mapping = Path('results/ancestral/baliphy-input-equivalence-20260927-v2/configuration_mapping.json')
    equivalence = mapping.parent / 'receipt.json'
    checked = json.loads(equivalence.read_text())
    assert sha(mapping) == checked['artifacts'][mapping.name]
    inputs = Path(baseline['input_receipt'])
    assert sha(inputs) == baseline['pins'][str(inputs)]
    original_jobs = {j['job_id']: j for j in json.loads(inputs.read_text())['jobs']}
    groups = collections.defaultdict(list)
    for row in json.loads(mapping.read_text()):
        groups[row['group']].append(row)
    assert len(groups) == 135 and sum(map(len, groups.values())) == 324
    jobs = []
    for group, rows in sorted(groups.items()):
        assert len({r['tree_sha256'] for r in rows}) == 1
        assert len({r['ordered_sequence_sha256'] for r in rows}) == 1
        representative = original_jobs[rows[0]['job_id']]
        for row in rows:
            job = original_jobs[row['job_id']]
            assert sha(job['alignment']) == row['alignment_sha256']
            assert sha(job['tree']) == row['tree_sha256']
        for label, prior in [('package', 'LogLaplace(6,2)'),
                             ('centered', 'LogLaplace(0,1)'),
                             ('broad', 'LogLaplace(0,2)')]:
            jobs.append(dict(representative,
                job_id=group + '-' + label, effective_input_group=group,
                original_configuration_ids=[r['job_id'] for r in rows],
                prior_label=label, alpha_prior=prior,
                smodel='LG +> F +> ASRV.Gamma(4,alpha=~' + prior + ')',
                imodel='RS07(rate=~LogLaplace(-4,0.707),meanLength=~ShiftedExponential(10,1))'))
    pins = dict(baseline['pins'])
    paths = [original, mapping, equivalence, Path(__file__),
             Path('scripts/run_baliphy_prior_initialization.py'),
             Path('metadata/baliphy_prior_definition_audit_20260927.json')]
    pins.update({str(p): sha(p) for p in paths})
    plan = dict(binary=baseline['binary'], jobs=jobs,
        output='results/ancestral/baliphy-prior-initialization-20260927-v1',
        options=['--test', '--seed', '20260929', '--indel-rates', 'constant'],
        timeout_seconds=600, memory_limit_bytes=6 * 1024**3,
        resources=dict(workers=2, cpus=2, memory_gib=16, swap_gib=0,
                       output_gib=5, planning_hours=[0.1, 34], paid_cost=0),
        pins=pins,
        scope='Initialization only: all135 distinct effective inputs, three alpha priors, '
              '405 executions with all324 original configuration aliases retained. '
              'Fixed trees and scale, free indel rate/length and amino-acid frequencies. '
              'No convergence or production posterior claim; aliases are not independent chains. '
              'Upper runtime estimate is 405 times 600 seconds divided by two workers, '
              'plus overhead. No GPU or paid resources.')
    target = Path('metadata/baliphy_prior_initialization_plan_20260927.json')
    assert not target.exists()
    target.write_text(json.dumps(plan, indent=2) + '\n')
    print(len(jobs), 'initializations prepared')


if __name__ == '__main__':
    main()
