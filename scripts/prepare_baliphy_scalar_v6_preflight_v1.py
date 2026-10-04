#!/usr/bin/env python3
"""Prepare the complete unchanged-model V6 native startup grid, without MCMC."""
import argparse
import ast
from datetime import datetime, timezone
import json
from pathlib import Path
import psutil

from ancestral_chain_attempt import sha, write_json
from baliphy_joint_sampler_scalar_v6 import build_jobs
from reference_measurement_union_sources import bind, verify


def project_sources(pins, entrypoints):
    pending = list(entrypoints)
    seen = set()
    while pending:
        path = pending.pop().resolve()
        if path in seen:
            continue
        seen.add(path)
        bind(pins, path)
        for node in ast.walk(ast.parse(path.read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom) and node.module
                     else [x.name for x in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                local = Path('scripts') / (name.split('.')[0] + '.py')
                if local.exists():
                    pending.append(local)
    return seen


def build(roles, previous_jobs, forbidden):
    samplers = build_jobs(roles, previous_jobs, forbidden)
    previous = {j['chain']['chain_id']: j for j in previous_jobs}
    jobs = []
    for sampler in samplers:
        future = sampler['chain']
        old = previous[sampler['source_chain_id']]['chain']
        command = sampler['config']['command']
        binary = command[command.index('--') + 1]
        config = dict(command=['/usr/bin/prlimit', '--as=' + str(12 * 2**30),
            '--cpu=600', '--fsize=' + str(256 * 2**20), '--', binary,
            '--seed', str(future['seed']), 'run', future['program'],
            '--test', '--log-format', 'json'], timeout_seconds=900,
            pins=dict(sampler['config']['pins']))
        chain = dict(future, seed=old['seed'], program=old['program'],
                     program_sha256=old['program_sha256'])
        jobs.append(dict(chain=chain, fresh_seed=future['seed'],
            source_chain_id=old['chain_id'], program=future['program'],
            project_scalar_schema=sampler['project_scalar_schema'], config=config))
    assert len(jobs) == 1620
    assert len({j['fresh_seed'] for j in jobs}) == 1620
    assert all(j['config']['command'][-3:] == ['--test', '--log-format', 'json']
               and '--iterations' not in j['config']['command'] for j in jobs)
    return jobs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--validation', type=Path, required=True)
    p.add_argument('--transport', type=Path, required=True)
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args()
    assert not a.plan.exists()
    gate = json.loads(a.validation.read_text())
    transport = json.loads(a.transport.read_text())
    assert gate['status'] == 'passed_full_scalar_v6_startup_software_contracts_v1'
    assert gate['full_roles'] == 1620 and gate['native_startups'] == 3
    assert gate['synthetic_full_producer_reader_serialization_checked']
    assert gate['artificial_failures_retained'] == 2
    assert transport['status'] == 'verified_original_scalar_v6_startup_software_wait_zero'
    assert transport['validation_sha256'] == sha(a.validation)
    verify(transport['source_hashes'])
    future_path = Path('metadata/baliphy_scalar_json_v6_future_models_20261004_v1.json')
    previous_path = Path('metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json')
    future = json.loads(future_path.read_text())
    previous = json.loads(previous_path.read_text())
    verify(future['pins']); verify(previous['pins'])
    roles = json.loads(Path(future['future_roles']).read_text())
    originals = json.loads(Path(previous['jobs']).read_text())
    forbidden = set(gate['forbidden_seeds'])
    assert len(forbidden) == future['forbidden_seed_count']
    assert set(gate['native_fixture_seeds']).isdisjoint(r['chain']['seed'] for r in roles)
    jobs = build(roles, originals, forbidden)
    root = a.inputs.resolve(); root.mkdir(exist_ok=False)
    job_path = root / 'jobs.json'; write_json(job_path, jobs)
    pins = dict(transport['source_hashes'])
    for mapping in [future['pins'], previous['pins']]:
        for name, digest in mapping.items(): bind(pins, name, digest)
    for path in [a.validation, a.transport, future_path, previous_path,
                 Path(future['future_roles']), Path(previous['jobs']), job_path]: bind(pins, path)
    sources = project_sources(pins, [Path(__file__),
        Path('scripts/run_baliphy_scalar_v6_preflight_v1.py'),
        Path('scripts/readback_baliphy_scalar_v6_preflight_v1.py')])
    closed_path = Path('metadata/baliphy_joint_logger_preflight_v5_completed_20261003.json')
    closed = json.loads(closed_path.read_text())
    archive = Path(closed['full_hash_archive'])
    assert sha(archive) == closed['full_hash_archive_sha256']
    proof = json.loads(archive.read_text())
    history = Path(closed['producer_receipt']).parent / 'dispositions.json'
    assert sha(history) == proof['source_hashes'][str(history)]
    rows = json.loads(history.read_text()); assert len(rows) == 1620
    seconds = sum(r['elapsed_seconds'] for r in rows)
    for path in [closed_path, archive, history]: bind(pins, path)
    verify(pins)
    resources = dict(checked_utc=datetime.now(timezone.utc).isoformat(), workers=2,
        cpus=2, memory_gib=32, swap_gib=0, blas_threads=1,
        per_worker_address_space_gib=12, per_worker_cpu_seconds_cap=600,
        per_worker_wall_seconds_cap=900, per_file_limit_mib=256,
        minimum_available_ram_gib=32, minimum_free_disk_gib=64, output_allowance_gib=16,
        historical_v5_full_startup_worker_seconds=seconds,
        sensitivity_worker_seconds=[seconds / 2, seconds * 4],
        sum_wall_timeout_worker_hours=1620 * 900 / 3600,
        available_memory_gib=psutil.virtual_memory().available / 2**30,
        available_disk_gib=psutil.disk_usage('.').free / 2**30,
        runtime_uncalibrated=True, finish_eta=None, native_startup_only=True,
        posterior_sampling=False, gpu=False, new_cost_usd=0,
        caveat='Historical startup work is not an ETA. Two12GiB AS caps plus8GiB '
               'group headroom. Output allowance is planning; native per-file limits '
               'and actual free-disk guard are separate.')
    plan = dict(jobs=str(job_path),
        output='results/ancestral/full-baliphy-scalar-v6-preflight-20261004-v1',
        completion='metadata/baliphy_scalar_v6_preflight_completed_20261004_v1.json',
        resources=resources, pins=pins, transitive_project_source_modules=len(sources),
        dependencies=[],
        scope='All1620preparedfreshV6roles,405quartets,135inputs,324aliases. Native--test '
              'initialization only: exact scalar inverse, reference homology, fixed-tip/free-ancestor '
              'representation, degree-aware density and rooted-tree checks. Strict initial CJSON '
              'records are distinct from MCMC V6 scalar rows. Retain every invalid/nativefailure '
              'role. Full original producer/readback/source/artifact/two-journal closure required. '
              'No MCMC, joint/scalar sampling qualification, posterior acceptance, original '
              'restart, GPU or paidresources.')
    with a.plan.open('x') as f: json.dump(plan, f, indent=2); f.write('\n')
    print(json.dumps(dict(plan=str(a.plan), **resources), indent=2))


if __name__ == '__main__':
    main()
