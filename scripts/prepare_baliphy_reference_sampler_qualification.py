#!/usr/bin/env python3
"""Prepare the full short native grid from the frozen startup and resource census."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

import psutil

from ancestral_chain_attempt import sha, write_json
from baliphy_reference_sampler_qualification import build_jobs, ITERATIONS
from launch_full_triad_sequence_geometry_followup import create
from run_baliphy_reference_preflight import verify


SCOPE = ('Full135effectiveinputs x3originalpriors x4freshseedroles, all324aliases: '
    '20 native sampler iterations per role after all1620correctedstartup records and their '
    'full provenance closure pass. No biological subset or pilot. Starting-state geometry '
    'changes only; original probability distribution, kernels, priors, observed unaligned '
    'tips and fixed-tip/free-ancestor representation retained. Family-wide48GiBaddress-space '
    'caps for historical allocation-risk families and12GiBotherwise;192GiBFIFOreservations '
    'held through native run and output audit under16CPU/200GiBcgroup/no swap. Every failed '
    'outcome retained without automatic retry; no posterior mixing, model/tree/root acceptance '
    'or full site-property replay implied. Full saved alignments/tips/candidate-node mapping, '
    'scalar logs, receipt/configuration/process/artifact/source hashes, reservation accounting '
    'and two original completion journals required for computational qualification closure. '
    'Original startup seeds now used for first MCMC: future independent longer sampling needs '
    'a new seed namespace disjoint from both original chains and this qualification. '
    'No GPU or new charges. All eight aims remain incomplete.')


def inputs():
    startup_path = Path('metadata/baliphy_reference_preflight_plan_20261003.json')
    startup = json.loads(startup_path.read_text()); verify(startup)
    closed_path = Path('metadata/baliphy_horizon_resources_completed_20261003.json')
    closed = json.loads(closed_path.read_text())
    assert closed['status'] == 'complete_verified_full_baliphy_horizon_resource_inventory'
    archive = Path(closed['full_hash_archive']); assert sha(archive) == closed['full_hash_archive_sha256']
    proof = json.loads(archive.read_text()); assert len(proof['services']) == 2
    census_paths = [archive.parent / (name + '.jsonl') for name in ['chains', 'attempts']]
    for p in census_paths: assert sha(p) == proof['source_hashes'][str(p)]
    source = Path('data/software_audits/baliphy-4.3-reference-initialization-20261003/source_manifest.json')
    manifest = json.loads(source.read_text())
    apis = [Path(x['installed_path']).resolve() for x in manifest['files'].values()]
    for p in apis: assert str(p) in startup['pins'] and sha(p) == startup['pins'][str(p)]
    binary = Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/bali-phy').resolve()
    paths = [startup_path, closed_path, archive, *census_paths, source]
    jobs, risk = build_jobs(json.loads(Path(startup['jobs']).read_text()),
        *[[json.loads(x) for x in p.read_text().splitlines()] for p in census_paths],
        binary, Path('/usr/bin/prlimit'), apis)
    return startup, jobs, risk, paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validation', type=Path, required=True); args = parser.parse_args()
    gate = json.loads(args.validation.read_text())
    assert gate['status'] == 'passed_full_reference_sampler_qualification_software_contracts'
    assert gate['full_chain_configurations_checked'] == gate['stress_roles'] == 1620
    assert gate['native_fixture_priors'] == ['broad', 'centered', 'package']
    assert gate['native_saved_alignments_checked'] == 9 and gate['native_candidate_frames_checked'] == 36
    assert gate['admission_abort_blocks_waiters'] and gate['fifo_large_role_not_bypassed']
    assert gate['failed_header_only_trace_retained'] and gate['failed_malformed_trace_retained']
    verify(dict(pins=gate['source_hashes']))
    startup, jobs, risk, paths = inputs()
    output = Path('results/ancestral/baliphy-reference-sampler-inputs-20261003-v1').resolve()
    output.mkdir(exist_ok=False); write_json(output / 'jobs.json', jobs)
    resources = dict(checked_utc=datetime.now(timezone.utc).isoformat(), workers=16, cpus=16,
        memory_gib=200, reservation_capacity_gib=192, swap_gib=0, blas_threads=1,
        risk_families=risk, role_address_space_gib_counts=dict(Counter(str(j['memory_reservation_bytes']//2**30) for j in jobs)),
        per_file_limit_gib=2, output_allowance_gib=64, minimum_free_disk_gib=128,
        iterations=ITERATIONS, saved_alignments_per_successful_role=3, candidate_frames_per_successful_role=12,
        linear_old_successful_worker_seconds=sum(j['linear_twenty_iteration_worker_seconds'] or 0 for j in jobs),
        sum_wall_timeout_worker_hours=sum(j['config']['timeout_seconds'] for j in jobs)/3600,
        available_memory_gib=psutil.virtual_memory().available/2**30,
        available_disk_gib=shutil.disk_usage(output).free/2**30, runtime_uncalibrated=True, finish_eta=None,
        posterior_qualified=False, gpu=False, new_cost_usd=0,
        caveat='Historical1000iteration success only, scaled20/1000; new starting state/launch/output overhead and two failures not calibrated. Reservations are declared caps, not empirical peak memory. Planning output allowance is not a measured guarantee or global hard quota; per-file cap/free-disk guard enforced. Full-grid short qualification is not sufficient posterior sampling.')
    owned = ['reference_sampler_memory_budget', 'baliphy_reference_sampler_qualification',
        'run_baliphy_reference_sampler_qualification', 'prepare_baliphy_reference_sampler_qualification',
        'launch_baliphy_reference_sampler_qualification', 'check_baliphy_reference_sampler_qualification',
        'record_baliphy_reference_sampler_checkpoint']
    paths += [Path('scripts/'+name+'.py') for name in owned]
    paths += [args.validation, output/'jobs.json', Path('metadata/baliphy_reference_startup_footer_plan_20261003.json'),
        Path('scripts/readback_independent_baliphy_chain.py'), Path('scripts/audit_baliphy_sample_mapping.py'),
        Path('scripts/baliphy_horizon_resource_inventory.py'),
        Path('results/ancestral/case-local-trees-20260927-v1/ancestral_node_mapping.tsv')]
    pins = dict(startup['pins']); pins.update({str(p):sha(p) for p in paths})
    plan = dict(jobs=str(output/'jobs.json'), output='results/ancestral/full-baliphy-reference-sampler-qualification-20261003-v1',
        startup_plan='metadata/baliphy_reference_startup_footer_plan_20261003.json',
        mapping='results/ancestral/case-local-trees-20260927-v1/ancestral_node_mapping.tsv',
        completion='metadata/baliphy_reference_sampler_qualification_completed_20261003.json',
        completion_plan='metadata/baliphy_reference_sampler_qualification_completion_plan_20261003.json',
        launch_inventory='metadata/baliphy_reference_sampler_qualification_launches_20261003.json',
        dependencies=['metadata/baliphy_reference_startup_footer_closure_launch_20261003.json'],
        software_validation=str(args.validation), pins=pins, resources=resources, scope=SCOPE)
    path = 'metadata/baliphy_reference_sampler_qualification_plan_20261003.json'; create(path, plan)
    print(json.dumps(dict(plan=path, **resources)), flush=True)


if __name__ == '__main__': main()
