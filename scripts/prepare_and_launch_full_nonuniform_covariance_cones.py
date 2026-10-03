#!/usr/bin/env python3
"""Run all closed operator cones with a distinct, unspecified residual diagonal."""
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

import psutil

from ancestral_chain_attempt import sha
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from reference_measurement_union_sources import verify
from run_full_nonuniform_covariance_cones import PRODUCER, READER, SUMMARY


def main():
    proofs = {}
    for name, version, session, status in [
        ('nonuniform_covariance_cone', 'v1', 85863,
         'passed_positive_diagonal_covariance_cone_dense_likelihood_contracts_v1'),
        ('full_nonuniform_covariance_cone', 'v2', 85487,
         'passed_complete_declared_nonuniform_covariance_cone_source_export_contracts_v2'),
    ]:
        paths = ['metadata/' + name + '_software_' + kind + '_20261003_' + version + '.json'
                 for kind in ['validation', 'execution', 'transport']]
        gate, execution, transport = [json.loads(Path(p).read_text()) for p in paths]
        assert gate['status'] == status and gate['scientific_eligibility'] is False
        assert execution['exit_code'] == 0 and execution['status'] == 'exited_zero_with_receipt'
        assert execution['receipt_sha256'] == sha(paths[0]) and not execution['timed_out']
        assert transport['actual_tool_session_id'] == session and transport['actual_tool_terminal_exit_code'] == 0
        assert transport['invocation_id'] == execution['invocation_id']
        assert transport['wrapper'] == execution['wrapper']
        assert transport['exact_wrapper_pid_journal_entries'] == transport['original_start_records'] == 1
        assert transport['original_completion_records'] == 1
        for d in [gate, execution, transport]:
            verify(d['source_hashes']); verify(d.get('artifacts', {}))
        if name.startswith('full_'):
            assert (gate['synthetic_cohorts'], gate['synthetic_certificates']) == (5, 10)
            assert len(gate['rehashed_output_cases_rejected']) == 12
            assert len(gate['rehashed_source_cases_rejected']) == 6
            assert gate['byte_exact_source_and_positive_artifact_restoration'] is True
        else:
            assert (gate['synthetic_diagonal_cases'], gate['covariance_cone_roundtrips'],
                    gate['ml_reml_dense_spectral_permutation_cases'], gate['finite_difference_score_cells']) == (20, 610, 120, 360)
            assert gate['malformed_certificate_cases_rejected'] == 56 and gate['invalid_diagonal_cases_rejected'] == 20
        proofs.update({p: sha(p) for p in paths})
    parent_plan = 'metadata/full_exact_covariance_folds_plan_20261003_v2.json'
    parent_completion = 'metadata/full_exact_covariance_folds_completed_20261003_v2.json'
    c = json.loads(Path(parent_completion).read_text())
    assert c['status'] == 'complete_verified_full_exact_uniform_covariance_folds_v2'
    assert (c['logical_cases'], c['cohorts'], c['certificates']) == (75188, 4340, 8680)
    assert sha(c['full_hash_archive']) == c['full_hash_archive_sha256']
    assert c['exact_process_journals_checked'] == 2 and c['scientific_eligibility'] is False
    resources = dict(checked_utc=datetime.now(timezone.utc).isoformat(), cpus=2, memory_gib=16,
        swap_gib=0, blas_threads=1, address_space_gib=12, cpu_seconds_per_stage=3600,
        per_file_limit_mib=256, output_allowance_gib=1, minimum_free_disk_gib=64,
        original_certificates=8680, original_cohorts=4340, forward_dimensions=[5, 6],
        original_variance_components=9, maximum_fraction_identity_products=8680 * 6 * 6 * 9,
        per_certificate_record_budget_bytes=16384,
        output_record_budget_bytes=8680 * 16384,
        runtime_planning_seconds_per_stage=[1, 600], runtime_calibrated=False,
        available_ram_gib=psutil.virtual_memory().available / 2**30,
        free_disk_gib=shutil.disk_usage('.').free / 2**30, new_cost_usd=0, gpu=False,
        scope='Full exact operator-cone metadata proof, not weighted data generation or fitting. '
              'Planning range is uncalibrated; one worker/one BLAS thread, CPU/AS/per-file and cgroup '
              'caps enforced. Record allowance is not a total directory quota or measured native peak. '
              'Full closed parent hashes and every original certificate are checked independently.')
    assert psutil.virtual_memory().available >= 16 * 2**30
    assert shutil.disk_usage('.').free >= 64 * 2**30
    rp = 'metadata/full_nonuniform_covariance_cone_resources_20261003_v1.json'; create(rp, resources)
    modules = ['nonuniform_covariance_cone', 'independent_nonuniform_covariance_cone',
        'full_nonuniform_covariance_cone_sources', 'run_full_nonuniform_covariance_cones',
        'check_nonuniform_covariance_cone', 'check_full_nonuniform_covariance_cones_v2',
        'prepare_and_launch_full_nonuniform_covariance_cones', 'reduced_covariance_basis',
        'covariance_exact_folds_v2', 'launch_baliphy_reference_sampler_qualification',
        'run_after_verified_dependencies_v2', 'close_full_triad_sequence_stage',
        'record_completed_process_handoffs_v2']
    paths = ['scripts/' + name + '.py' for name in modules] + [rp, parent_plan, parent_completion,
        c['full_hash_archive'], '/usr/bin/prlimit', sys.executable, *proofs]
    output = 'results/phylogeny/full-positive-diagonal-covariance-cones-20261003-v1'
    assert not Path(output).exists()
    pp = 'metadata/full_nonuniform_covariance_cone_plan_20261003_v1.json'
    completion = 'metadata/full_nonuniform_covariance_cones_completed_20261003_v1.json'
    scope = ('Complete closed4340cohorts/8680mode certificates/75188logical cases/68220240row occurrences '
        'positive-diagonal covariance-cone proof. Keep residual D distinct from target I for any finite '
        'strictly positive diagonal; exact gene/model/family folds and genuine pair exceptions retained. '
        'Full original parent source hashes checked; no original operator/Gram/certificate edits. '
        'Independent Fraction column-image reconstruction and exact nonnegative right-inverse identities '
        'check every serialized record. Full source/artifact/two new original journals gate closure. '
        'No actual diagonal chosen, confidence-to-variance policy calibrated, uniform numerical envelope '
        'inherited, raw/REML weighted basis qualified, weighted fit, native restart, GPU, new charge or '
        'biological aim completed. Actual weighted inputs/identifiability/timing/controls remain separate.')
    create(pp, dict(parent_plan=parent_plan, parent_completion=parent_completion,
        expected=dict(logical_cases=75188, cohorts=4340, certificates=8680), output=output,
        resources=resources, pins={p: sha(p) for p in paths}, completion=completion, scope=scope))
    prefix = ['/usr/bin/prlimit', '--as=' + str(12 * 2**30), '--cpu=3600',
              '--fsize=' + str(256 * 2**20), '--']
    producer = launch('full-nonuniform-covariance-cones-v1',
        prefix + [sys.executable, 'scripts/run_full_nonuniform_covariance_cones.py', '--plan', pp],
        [], pp, cpus=2, memory=16)
    reader = launch('full-nonuniform-covariance-cones-v1-readback',
        prefix + [sys.executable, 'scripts/run_full_nonuniform_covariance_cones.py', '--plan', pp, '--reader'],
        [producer], pp, cpus=2, memory=16)
    cp = 'metadata/full_nonuniform_covariance_cone_completion_plan_20261003_v1.json'
    create(cp, dict(source_plan=pp, producer_receipt=output + '/receipt.json',
        independent_readback=output + '/readback.json', producer_status=PRODUCER, reader_status=READER,
        completed_status='complete_verified_full_positive_diagonal_covariance_cones_v1',
        summary_fields=SUMMARY, launches=[producer, reader],
        pins={p: sha(p) for p in [pp, producer, reader, 'scripts/close_full_triad_sequence_stage.py',
                                'scripts/record_completed_process_handoffs_v2.py']},
        output=completion, scope=scope))
    closer = launch('full-nonuniform-covariance-cones-v1-closure',
        [sys.executable, 'scripts/close_full_triad_sequence_stage.py', '--plan', cp],
        [producer, reader], cp, cpus=2, memory=16)
    create('metadata/full_nonuniform_covariance_cone_launches_20261003_v1.json',
        dict(status='queued_complete_closed_operator_positive_diagonal_cone_proof',
             source_plan=pp, source_plan_sha256=sha(pp), launches=[producer, reader, closer],
             expected_certificates=8680, residual_diagonal_prepared=False,
             weighted_qualification_complete=False, production_fitting_launched=False,
             gpu=False, new_cost_usd=0, all_eight_aims_incomplete=True))


if __name__ == '__main__':
    main()
