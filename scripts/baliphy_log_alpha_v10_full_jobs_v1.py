"""Construct all 24 V10 roles from the closed original V7 comparison matrix."""
import copy
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_joint_fasta_v7_failure_jobs_v1 import (
    SUFFIX as V7_SUFFIX, failed_sources, validate_jobs as validate_v7_jobs,
)
from baliphy_log_alpha_logger_v10 import transform, reverse
from reference_measurement_union_sources import bind, verify


SUFFIX = '-latent-log-alpha-v10-full-comparison'


def sources():
    original, pins, mapping = failed_sources()
    plan_path = Path('metadata/baliphy_joint_fasta_v7_failure_grid_plan_20261004_v1.json')
    completed_path = Path('metadata/baliphy_joint_fasta_v7_failure_grid_completed_20261005_v1.json')
    plan, completed = [json.loads(p.read_text()) for p in (plan_path, completed_path)]
    assert completed['status'] == 'complete_verified_all24_v7_failure_comparison_dispositions'
    logger_path = Path('metadata/baliphy_log_alpha_v10_completed_20261005_v1.json')
    logger = json.loads(logger_path.read_text())
    assert logger['status'] == 'complete_verified_V10_latent_alpha_logger_software_and_diagnostic_replay'
    assert logger['full_generated_sources_reversibly_checked'] == 405
    assert logger['paired_native_prior_controls'] == 3
    assert logger['unchanged_original_scientific_files'] == 18
    jobs_path = Path(plan['jobs'])
    old = json.loads(jobs_path.read_text())
    scope = validate_v7_jobs(old, original)
    dispositions_path = Path(plan['output']) / 'dispositions.json'
    dispositions = json.loads(dispositions_path.read_text())
    by_id = {r['chain_id']: r for r in dispositions}
    assert len(dispositions) == len(by_id) == 24
    assert set(by_id) == {j['chain']['chain_id'] for j in old}
    assert all(r['exit_code'] == 0 for r in dispositions)
    for job in old:
        row = by_id[job['chain']['chain_id']]
        assert row['seed'] == job['chain']['seed']
        job['paired_v7_native_receipt'] = row['native_receipt']
        job['paired_v7_native_receipt_sha256'] = row['native_receipt_sha256']
        bind(pins, row['native_receipt'], row['native_receipt_sha256'])
    for mapping_ in [plan['pins'], completed['source_hashes'], logger['source_hashes']]:
        for path, digest in mapping_.items():
            bind(pins, path, digest)
    for path in [plan_path, completed_path, logger_path, jobs_path, dispositions_path]:
        bind(pins, path)
    verify(pins)
    return old, pins, mapping, scope


def candidate_job(previous, program):
    program = Path(program).resolve()
    old = previous['chain']
    assert old['chain_id'].endswith(V7_SUFFIX)
    assert sha(old['program']) == old['program_sha256']
    assert reverse(program.read_text()) == Path(old['program']).read_text()
    job = copy.deepcopy(previous)
    job['chain'].update(chain_id=old['chain_id'][:-len(V7_SUFFIX)] + SUFFIX,
                        program=str(program), program_sha256=sha(program))
    command = job['config']['command']
    assert command[command.index('run') + 1] == old['program']
    command[command.index('run') + 1] = str(program)
    assert job['config']['pins'].pop(old['program']) == old['program_sha256']
    job['config']['pins'][str(program)] = sha(program)
    job['source_v7_chain_id'] = old['chain_id']
    job['source_v7_program_sha256'] = old['program_sha256']
    return job


def build_jobs(previous, root):
    root = Path(root).resolve()
    root.mkdir(exist_ok=False)
    jobs, models = [], {}
    for job in previous:
        old = job['chain']
        key = old['effective_input_group'] + '-' + old['prior_label']
        program = root / (key + '.hs')
        text = transform(Path(old['program']).read_text())
        if key not in models:
            program.write_text(text)
            models[key] = program
        else:
            assert program.read_text() == text
        jobs.append(candidate_job(job, program))
    assert len(models) == 6
    return jobs


def validate_jobs(jobs, previous, original_scope):
    by_id = {j['chain']['chain_id']: j for j in previous}
    assert len(jobs) == len(previous) == len(by_id) == 24
    assert {j['source_v7_chain_id'] for j in jobs} == set(by_id)
    assert len({j['chain']['chain_id'] for j in jobs}) == 24
    assert len({j['chain']['seed'] for j in jobs}) == 24
    for job in jobs:
        assert job == candidate_job(by_id[job['source_v7_chain_id']], job['chain']['program'])
        verify(job['config']['pins'])
    assert len({j['chain']['program'] for j in jobs}) == 6
    return dict(original_scope, paired_v7_roles=24, v10_reversible_models=6,
                only_native_program_and_chain_namespace_changed=True,
                full_input_tips=622, full_V10_matrix_roles=24,
                scientific_eligibility=False, posterior_qualified=False)
