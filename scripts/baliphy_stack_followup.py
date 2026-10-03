"""Scoped stack correction, fresh seeds and complete original-role accounting."""
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_joint_sampler_gates_v3 import closed_stage, resource_boundaries
from baliphy_joint_sampler_qualification_v3 import SUCCESS, summarize
from reference_measurement_union_sources import verify


NAMESPACE = 'fungal-baliphy-stack-followup-20261003-v1'
DIAGNOSIS_STATUS = 'validated_two_input_stack_exhaustion_and_initial_joint_frames'
SUMMARY_FIELDS = ['full_chains','full_quartets','effective_inputs','original_configuration_aliases',
    'original_status_counts','original_failed_roles','followup_roles','followup_status_counts',
    'selected_checked_sampler_attempts','selected_unsuccessful_sampler_attempts','selected_complete_quartets',
    'selected_unresolved_quartets','selected_joint_saved_frames','selected_joint_ancestral_residue_category_pairs',
    'all_original_failures_retained','all_fresh_seeds_disjoint','posterior_qualified','resource_audit','reservation_audit']


def forbidden_seeds(originals, historical):
    seeds = {20266001,20266002,20266003,20267001,20267002,20267003}
    for job in originals + historical:
        seeds.update([job['chain']['seed'],job['source_seed']])
    return seeds


def build_jobs(originals, rows, historical, diagnosis):
    assert diagnosis['status'] == DIAGNOSIS_STATUS
    assert diagnosis['full_original_roles'] == 1620 and diagnosis['original_failed_roles'] == 24
    assert diagnosis['full_twenty_iteration_correction_qualified'] is diagnosis['posterior_qualified'] is False
    assert len(originals) == len(rows) == len(historical) == 1620
    by = {j['chain']['chain_id']:j for j in originals}; rb = {r['chain_id']:r for r in rows}
    assert len(by) == len(rb) == 1620 and set(by) == set(rb)
    assert len({j['chain']['effective_input_group'] for j in originals}) == 135
    assert len({a for j in originals for a in j['chain']['original_configuration_ids']}) == 324
    assert Counter(j['chain']['prior_label'] for j in originals) == {'broad':540,'centered':540,'package':540}
    groups = Counter((j['chain']['effective_input_group'],j['chain']['prior_label']) for j in originals)
    assert len(groups) == 405 and set(groups.values()) == {4}
    failed = [r for r in rows if r['status'] != SUCCESS]
    assert len(failed) == 24 and Counter(r['effective_input_group'] for r in failed) == diagnosis['failed_input_group_counts']
    blocked = forbidden_seeds(originals,historical); new = set(); jobs = []
    for row in sorted(failed,key=lambda r:r['chain_id']):
        old = by[row['chain_id']]; c = old['chain']
        assert row['exit_code'] == -11 and row['joint_frames'] == [] and row['saved_alignments'] == 0
        assert row['bad_alloc'] is False and row['allocation_warning_lines'] == 0
        assert c['family'] == 'OG0000972' and c['proteins'] == 622 and old['memory_reservation_bytes'] == 48*2**30
        for key,value in [('effective_input_group',c['effective_input_group']),('prior_label',c['prior_label']),
                          ('chain_role',c['chain']),('seed',c['seed']),('original_configuration_ids',c['original_configuration_ids'])]:
            assert row[key] == value
        ordinal = 0
        while True:
            token = (NAMESPACE + ':' + c['chain_id'] + ':' + str(ordinal)).encode()
            seed = 1 + int.from_bytes(hashlib.sha256(token).digest()[:8],'big') % (2**31-1)
            if seed not in blocked | new: break
            ordinal += 1
        new.add(seed)
        chain = dict(c,chain_id=c['chain_id']+'-stack64-v1',seed=seed)
        config = copy.deepcopy(old['config']); command = config['command']
        assert '--stack=67108864' not in command and not any(x.startswith('--stack=') for x in command)
        assert command[command.index('--seed')+1] == str(c['seed'])
        assert command[command.index('--iterations')+1] == '20' and command[command.index('run')+1] == c['program']
        command[command.index('--seed')+1] = str(seed); command.insert(command.index('--'),'--stack=67108864')
        # Removing the scoped stack flag and restoring the seed exactly
        # recovers the original command; every config pin and timeout survives.
        restored = [x for x in command if x != '--stack=67108864']
        restored[restored.index('--seed')+1] = str(c['seed']); assert restored == old['config']['command']
        assert config['pins'] == old['config']['pins'] and config['timeout_seconds'] == old['config']['timeout_seconds']
        jobs.append(dict(old,chain=chain,config=config,source_seed=c['seed'],source_chain_id=c['chain_id'],
            initialization_only_seed_reused_for_first_mcmc=False,original_failure_receipt=row['native_receipt'],
            original_failure_receipt_sha256=row['native_receipt_sha256'],seed_namespace=NAMESPACE,seed_collision_ordinal=ordinal))
    assert len(new) == len(jobs) == 24 and new.isdisjoint(blocked)
    assert Counter(j['chain']['prior_label'] for j in jobs) == {'broad':8,'centered':8,'package':8}
    assert len({(j['chain']['effective_input_group'],j['chain']['prior_label']) for j in jobs}) == 6
    return jobs


def prerequisites(plan):
    source,closed,bindings = closed_stage(plan['source_plan'],'complete_verified_full_joint_short_sampler_qualification_v3')
    assert (closed['full_chains'],closed['checked_sampler_attempts'],closed['unsuccessful_sampler_attempts']) == (1620,1596,24)
    assert (closed['full_quartets'],closed['effective_inputs'],closed['original_configuration_aliases']) == (405,135,324)
    observer,oc,ob = closed_stage(plan['source_resource_plan'],'complete_verified_full_joint_sampler_resource_observation_v3')
    assert observer['sampler_plan'] == plan['source_plan']; resource_boundaries(observer,oc,ob)
    for name,h in ob.items():
        assert name not in bindings or bindings[name] == h; bindings[name] = h
    diagnosis = json.loads(Path(plan['diagnosis']).read_text()); assert diagnosis['status'] == DIAGNOSIS_STATUS
    verify(diagnosis['source_hashes']); bindings[str(plan['diagnosis'])] = sha(plan['diagnosis'])
    originals = json.loads(Path(source['jobs']).read_text())
    rows_path = Path(source['output'])/'dispositions.json'
    assert sha(rows_path) == bindings[str(rows_path)]
    rows = json.loads(rows_path.read_text())
    historical_plan = json.loads(Path(source['historical_sampler_plan']).read_text())
    historical = json.loads(Path(historical_plan['jobs']).read_text())
    jobs = json.loads(Path(plan['jobs']).read_text())
    assert jobs == build_jobs(originals,rows,historical,diagnosis)
    assert plan['resources']['iterations'] == 20 and plan['resources']['native_stack_bytes'] == 64*2**20
    assert plan['resources']['workers'] == 4 and plan['resources']['reservation_capacity_gib'] == 192
    # Source closure already independently reconstructed all 1620 original
    # outputs. Rehash its complete archive; never regenerate original arrays.
    verify(bindings)
    return originals,rows,jobs,bindings


def merge(originals,rows,followups,followup_rows,source_root,root):
    assert len(rows) == len(originals) == 1620
    new = {r['chain_id']:r for r in followup_rows}; assert len(new) == len(followups) == 24
    follow = {j['source_chain_id']:j for j in followups}; assert len(follow) == 24
    selected = []; ledger = []
    for original in sorted(rows,key=lambda r:r['chain_id']):
        cid = original['chain_id']; source = Path(source_root)/'chains'/(cid+'.json')
        assert json.loads(source.read_text()) == original
        entry = dict(original_chain_id=cid,original_record=str(source),original_record_sha256=sha(source),
            original_status=original['status'],original_native_receipt=original['native_receipt'],
            original_native_receipt_sha256=original['native_receipt_sha256'],followup_record=None,
            followup_record_sha256=None,followup_status=None)
        if cid in follow:
            job = follow[cid]; fresh = new[job['chain']['chain_id']]
            for key in ['effective_input_group','prior_label','chain_role','family','original_configuration_ids']:
                assert fresh[key] == original[key]
            assert fresh['seed'] == job['chain']['seed'] and fresh['source_seed'] == original['seed']
            fp = Path(root)/'chains'/(fresh['chain_id']+'.json')
            assert json.loads(fp.read_text()) == fresh
            entry.update(followup_record=str(fp),followup_record_sha256=sha(fp),followup_status=fresh['status'])
            selected.append(dict(fresh,chain_id=cid,selected_followup_chain_id=fresh['chain_id']))
        else:
            assert original['status'] == SUCCESS
            selected.append(original)
        ledger.append(entry)
    science = summarize(selected)
    summary = dict(full_chains=1620,full_quartets=405,effective_inputs=135,original_configuration_aliases=324,
        original_status_counts=dict(Counter(r['status'] for r in rows)),original_failed_roles=len(follow),followup_roles=len(new),
        followup_status_counts=dict(Counter(r['status'] for r in new.values())),
        selected_checked_sampler_attempts=science['checked_sampler_attempts'],
        selected_unsuccessful_sampler_attempts=science['unsuccessful_sampler_attempts'],
        selected_complete_quartets=science['complete_quartets'],selected_unresolved_quartets=science['unresolved_quartets'],
        selected_joint_saved_frames=science['joint_saved_frames'],
        selected_joint_ancestral_residue_category_pairs=science['joint_ancestral_residue_category_pairs'],
        all_original_failures_retained=True,all_fresh_seeds_disjoint=True,posterior_qualified=False)
    assert len(ledger) == 1620 and sum(r['followup_record'] is not None for r in ledger) == 24
    return ledger,summary
