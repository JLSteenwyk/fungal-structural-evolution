"""Closed full-grid prerequisites for the future joint logging qualification."""
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from run_baliphy_reference_preflight import verify


def closed_stage(plan_path, status):
    path = Path(plan_path); plan = json.loads(path.read_text()); verify(plan)
    completed_path = Path(plan['completion']); completed = json.loads(completed_path.read_text())
    assert completed['status'] == status and completed['scientific_eligibility'] is False
    archive_path = Path(completed['full_hash_archive'])
    assert sha(archive_path) == completed['full_hash_archive_sha256']
    proof = json.loads(archive_path.read_text())
    assert len(proof['services']) == completed['exact_process_journals_checked'] == 2
    for service in proof['services']:
        assert service['captured_process_messages'] > 0 and service['completion_resource_records']
        assert service['observed_state']['Result'] == 'success' and service['observed_state']['ExecMainStatus'] == '0'
    bindings = dict(proof['source_hashes'])
    assert str(path) in bindings and bindings[str(path)] == sha(path)
    assert sha(completed['producer_receipt']) == completed['producer_receipt_sha256']
    assert sha(completed['independent_readback']) == completed['independent_readback_sha256']
    for name in ['producer_receipt','independent_readback']:
        receipt = json.loads(Path(completed[name]).read_text())
        assert receipt['plan_sha256'] == sha(path) and receipt['scientific_eligibility'] is False
    for p in [path,completed_path,archive_path]:
        assert str(p) not in bindings or bindings[str(p)] == sha(p)
        bindings[str(p)] = sha(p)
    return plan,completed,bindings


def resource_boundaries(plan, completed, bindings):
    assert completed['expected_roles'] == completed['roles_with_native_attempt'] == 1620
    assert completed['roles_without_native_attempt'] == 0
    assert completed['roles_with_live_observation'] + completed['roles_without_live_observation'] == 1620
    assert completed['sampler_controller_terminal_status'] == 'verified_original_terminal_success_with_bound_completed_artifacts'
    assert completed['posterior_qualified'] is False
    observed = completed['maximum_observed_cgroup_reported_peak_bytes']
    assert type(observed) is int and 0 <= observed <= 200*2**30
    path = Path(plan['output'])/'observations.jsonl'
    assert sha(path) == bindings[str(path)]
    count = 0
    with path.open() as handle:
        for count,line in enumerate(handle,1):
            snapshot = json.loads(line); assert snapshot['sequence'] == count-1
            group = snapshot['cgroup']
            assert group['limits'] == {'cpu.max':'1600000 100000','memory.max':str(200*2**30),'memory.swap.max':'0'}
            assert group['memory_bytes']['memory.swap.current'] == 0
            assert all(group['memory_events'][key] == 0 for key in ['max','oom','oom_kill','oom_group_kill'])
    assert count > 0
    # Sample gaps and roles without live observations remain in the closed
    # record. The observed counters are not final per-native memory bounds.


def prerequisites(plan):
    bindings = {}
    def merge(additional):
        for p,h in additional.items():
            assert p not in bindings or bindings[p] == h,p
            bindings[p] = h
    startup,closed,proof = closed_stage(plan['startup_plan'],'complete_verified_full_joint_logger_startup')
    merge(proof)
    assert (closed['full_chains'],closed['full_quartets'],closed['effective_inputs'],
        closed['original_configuration_aliases'],closed['validated_startups'],closed['unsuccessful_startups'],
        closed['complete_startup_quartets'],closed['unresolved_startup_quartets']) == (1620,405,135,324,1620,0,405,0)
    assert closed['posterior_sampling_launched'] is False
    path = Path(startup['output'])/'dispositions.json'; assert sha(path) == bindings[str(path)]
    rows = json.loads(path.read_text()); jobs = json.loads(Path(plan['jobs']).read_text())
    by_id = {j['chain']['chain_id']:j for j in jobs}
    assert len(rows) == len({r['chain_id'] for r in rows}) == len(by_id) == 1620
    assert set(by_id) == {r['chain_id'] for r in rows}
    for row in rows:
        job = by_id[row['chain_id']]; chain = job['chain']
        assert row['status'] == 'reference_startup_homology_density_and_representation_checked'
        assert row['fresh_seed'] == chain['seed'] and row['source_seed'] == job['source_seed']
        assert row['scientific_eligibility'] is row['posterior_sampling_launched'] is False
        for key,value in [('chain_role',chain['chain']),('prior_label',chain['prior_label']),
            ('effective_input_group',chain['effective_input_group']),('original_configuration_ids',chain['original_configuration_ids'])]:
            assert row[key] == value
    native,native_closed,native_proof = closed_stage(plan['historical_sampler_plan'],
        'complete_verified_full_reference_short_sampler_qualification'); merge(native_proof)
    assert native_closed['full_chains'] == native_closed['checked_sampler_attempts'] + native_closed['unsuccessful_sampler_attempts'] == 1620
    assert native_closed['full_quartets'] == 405 and native_closed['effective_inputs'] == 135
    assert native_closed['original_configuration_aliases'] == 324 and native_closed['posterior_qualified'] is False
    observer,observer_closed,observer_proof = closed_stage(plan['historical_resource_plan'],
        'complete_verified_full_sampler_resource_observation'); merge(observer_proof)
    assert observer['sampler_plan'] == plan['historical_sampler_plan']
    resource_boundaries(observer,observer_closed,observer_proof)
    reader,reader_closed,reader_proof = closed_stage(plan['historical_replay_plan'],
        'complete_verified_full_independent_short_sampler_native_replay_v2'); merge(reader_proof)
    assert reader['sampler_plan'] == plan['historical_sampler_plan']
    assert reader_closed['full_chains'] == reader_closed['checked_chains'] + reader_closed['unresolved_chains'] == 1620
    assert reader_closed['full_groups'] == 405 and reader_closed['effective_inputs'] == 135
    assert reader_closed['original_configuration_aliases'] == 324
    assert reader_closed['posterior_qualified'] is reader_closed['ancestral_categories_available'] is reader_closed['tip_logs_joint_trajectory_available'] is False
    verify({'pins':bindings})
    return bindings
