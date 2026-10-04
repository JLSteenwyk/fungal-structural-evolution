"""Closed V6 startup plus the preserved original historical prerequisite checks."""
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_joint_sampler_gates_v3 import closed_stage, prerequisites as historical_prerequisites
from baliphy_scalar_json_logger_v6c import SCHEMA
from reference_measurement_union_sources import bind, verify


def startup_roles(startup, closed, bindings, jobs):
    assert (closed['full_chains'], closed['full_quartets'], closed['effective_inputs'],
        closed['original_configuration_aliases'], closed['validated_startups'],
        closed['unsuccessful_startups'], closed['complete_startup_quartets'],
        closed['unresolved_startup_quartets']) == (1620, 405, 135, 324, 1620, 0, 405, 0)
    assert closed['posterior_sampling_launched'] is False
    path = Path(startup['output']) / 'dispositions.json'
    assert sha(path) == bindings[str(path)]
    rows = json.loads(path.read_text())
    by_id = {j['chain']['chain_id']: j for j in jobs}
    assert len(rows) == len({r['chain_id'] for r in rows}) == len(by_id) == len(jobs) == 1620
    assert {r['chain_id'] for r in rows} == set(by_id)
    for row in rows:
        job = by_id[row['chain_id']]; chain = job['chain']
        assert row['status'] == 'reference_startup_homology_density_and_representation_checked'
        assert row['project_scalar_schema'] == job['project_scalar_schema'] == SCHEMA
        assert row['fresh_seed'] == chain['seed'] and row['source_seed'] == job['source_seed']
        assert row['scientific_eligibility'] is row['posterior_sampling_launched'] is False
        for key, value in [('chain_role', chain['chain']), ('prior_label', chain['prior_label']),
            ('effective_input_group', chain['effective_input_group']),
            ('original_configuration_ids', chain['original_configuration_ids'])]:
            assert row[key] == value


def prerequisites(plan):
    # Missing/unfinished startup refuses before any new native role is admitted.
    startup, closed, bindings = closed_stage(plan['startup_plan'],
        'complete_verified_full_scalar_v6_startup_v1')
    jobs = json.loads(Path(plan['jobs']).read_text())
    startup_roles(startup, closed, bindings, jobs)
    previous_path = Path(plan['previous_joint_sampler_plan'])
    previous = json.loads(previous_path.read_text()); verify(previous['pins'])
    # Keep all four original V3 prerequisites: historical native/source custody,
    # memory observation, independent sequence-coordinate/rate-review accounting,
    # and the earlier startup closure. These never accept historical numbers.
    for name, digest in historical_prerequisites(previous).items(): bind(bindings, name, digest)
    old_jobs = json.loads(Path(previous['jobs']).read_text())
    old_by_id = {j['chain']['chain_id']: j for j in old_jobs}
    assert len(old_jobs) == len(old_by_id) == 1620
    assert {j['source_chain_id'] for j in jobs} == set(old_by_id)
    for job in jobs:
        old = old_by_id[job['source_chain_id']]['chain']; chain = job['chain']
        assert job['source_seed'] == old['seed']
        for field in ['effective_input_group', 'prior_label', 'chain', 'original_configuration_ids',
                      'family', 'proteins', 'alignment', 'alignment_sha256', 'tree', 'tree_sha256']:
            assert chain[field] == old[field]
    for path in [previous_path, Path(previous['jobs']), Path(plan['jobs'])]: bind(bindings, path)
    verify(bindings)
    return bindings
