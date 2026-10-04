#!/usr/bin/env python3
"""Prepare the exact qualified V6 computational grid and resource/closure plans."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

import psutil

from ancestral_chain_attempt import sha, write_json
from baliphy_scalar_json_logger_v6c import SCHEMA
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind, verify
from run_baliphy_scalar_v6_resource_observer_v1 import qualified_jobs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.plan.exists() and not args.inputs.exists() and not args.receipt.exists()
    pins = {}
    transports = {
        'sampler_software_transport':'metadata/baliphy_scalar_v6_sampler_controller_software_transport_20261004_v1.json',
        'observer_software_transport':'metadata/baliphy_scalar_v6_resource_observer_software_transport_20261004_v1.json'}
    expected = ['verified_original_scalar_v6_sampler_controller_software_wait_zero',
                'verified_original_scalar_v6_resource_observer_software_wait_zero']
    for (field, name), status in zip(transports.items(), expected):
        value = json.loads(Path(name).read_text())
        assert value['status'] == status and value['actual_tool_terminal_exit_code'] == 0
        assert value['entire_terminal_payload_matched']
        verify(value['source_hashes'])
        for path, digest in value['source_hashes'].items(): bind(pins, path, digest)
        bind(pins, name)
    future_path = Path('metadata/baliphy_scalar_json_v6_future_models_20261004_v1.json')
    previous_path = Path('metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json')
    startup_path = Path('metadata/baliphy_scalar_v6_preflight_plan_20261004_v1.json')
    inventory_path = Path('metadata/baliphy_scalar_v6_preflight_launches_20261004_v1.json')
    future = json.loads(future_path.read_text()); previous = json.loads(previous_path.read_text())
    inventory = json.loads(inventory_path.read_text())
    assert inventory['source_plan_sha256'] == sha(startup_path)
    for mapping in [future['pins'], previous['pins']]:
        verify(mapping)
        for path, digest in mapping.items(): bind(pins, path, digest)
    reference = Path('data/software_audits/baliphy-scalar-json-v6-sampler-controller-20261004-v1/actual_full_grid_jobs.json').resolve()
    jobs = qualified_jobs(dict(jobs=str(reference), scalar_schema=SCHEMA))
    roles = json.loads(Path(future['future_roles']).read_text())
    by_id = {row['chain']['chain_id']:row for row in roles}
    assert len(by_id) == len(jobs) == 1620
    assert all(job['chain'] == by_id[job['chain']['chain_id']]['chain'] for job in jobs)
    assert Counter(job['memory_reservation_bytes']//2**30 for job in jobs) == {12:1512,48:108}
    assert len({job['chain']['effective_input_group'] for job in jobs}) == 135
    assert len({alias for job in jobs for alias in job['chain']['original_configuration_ids']}) == 324
    assert len({job['chain']['seed'] for job in jobs}) == 1620
    # This copies the qualified definitions, not a fresh census of artificial
    # telemetry configurations. Artificial fixture seeds are not native usage.
    args.inputs.mkdir(parents=True, exist_ok=False)
    job_path = (args.inputs/'jobs.json').resolve()
    with job_path.open('xb') as handle: handle.write(reference.read_bytes())
    assert sha(job_path) == sha(reference)
    dependencies = [inventory['launches'][2], *previous['dependencies']]
    for path in [future_path, previous_path, startup_path, inventory_path, reference, job_path,
                 Path(future['future_roles']), Path(previous['mapping']), *map(Path, dependencies)]: bind(pins, path)
    modules = project_sources(pins, [Path(__file__), Path('scripts/launch_baliphy_scalar_v6_sampler_execution_v1.py'),
        Path('scripts/run_baliphy_scalar_v6_sampler_v1.py'), Path('scripts/run_baliphy_scalar_v6_resource_observer_v1.py')])
    resources = dict(previous['resources'])
    resources.update(checked_utc=datetime.now(timezone.utc).isoformat(),
        available_memory_gib=psutil.virtual_memory().available/2**30,
        available_disk_gib=psutil.disk_usage('.').free/2**30,
        runtime_uncalibrated=True, finish_eta=None,
        caveat='Full20iteration computational qualification of all1620V6roles, not adequate posterior. '
               'Original per-role timeouts/12GiB1512roles/48GiB108roles retained. Sixteen workers/192GiB '
               'declared AS leases under200GiB/noSwap with8GiB group headroom; these are not observed '
               'final native peaks or a long-chain bound. Historical worker cost is not an ETA. '
               '128GiB output planning and256GiB free-disk guard; no hard global output quota. '
               'Observer is armed concurrently; actual telemetry/source/artifact closure follows native termination.')
    scope = ('All1620qualified V6 roles/405quartets/135inputs/324aliases: unchanged20iteration '
             'computational gate with preserved model/priors/initialization and exact future seeds. '
             'Full original V6startup plus four historical source/memory/sequence-coordinate-review '
             'prerequisites close before native admission. Historical numeric values never adopted. '
             'Malformed/nonfinite/native-failure roles and unresolved quartets remain explicit, no '
             'automatic retry. Original independent serialized replay and full source/artifact/two '
             'original-journal closure mandatory. Concurrent read-only observer keeps missing readings '
             'and requires its own closure after native termination. No adequate posterior, global '
             'model/root/predictor validation, biological acceptance, original restart, GPU or charges.')
    plan = dict(jobs=str(job_path), scalar_schema=SCHEMA,
        output='results/ancestral/full-baliphy-scalar-v6-short-sampler-20261004-v1', mapping=previous['mapping'],
        startup_plan=str(startup_path), startup_closure_launch=inventory['launches'][2],
        previous_joint_sampler_plan=str(previous_path), dependencies=dependencies,
        completion='metadata/baliphy_scalar_v6_short_sampler_completed_20261004_v1.json',
        completion_plan='metadata/baliphy_scalar_v6_short_sampler_completion_plan_20261004_v1.json',
        launch_inventory='metadata/baliphy_scalar_v6_sampler_execution_launches_20261004_v1.json',
        observer_output='results/ancestral/full-baliphy-scalar-v6-resource-observation-20261004-v1',
        observer_plan='metadata/baliphy_scalar_v6_resource_observation_plan_20261004_v1.json',
        observer_completion='metadata/baliphy_scalar_v6_resource_observation_completed_20261004_v1.json',
        observer_completion_plan='metadata/baliphy_scalar_v6_resource_observation_completion_plan_20261004_v1.json',
        **transports, pins=pins, resources=resources, transitive_project_source_modules=len(modules), scope=scope)
    assert not Path(plan['output']).exists() and not Path(plan['observer_output']).exists()
    verify(pins); write_json(args.plan, plan)
    bind(pins, args.plan)
    result = dict(status='prepared_exact_qualified_full_scalar_v6_sampler_execution',
        checked_utc=datetime.now(timezone.utc).isoformat(), plan=str(args.plan), plan_sha256=sha(args.plan),
        native_roles=1620, quartets=405, effective_inputs=135, aliases=324,
        copied_job_file_matches_qualified_reference=True, expected_original_prerequisite_closures=5,
        actual_startup_closure_present=Path(json.loads(startup_path.read_text())['completion']).exists(),
        source_hashes=pins, transitive_project_source_modules=len(modules),
        native_execution_launched=False, posterior_qualified=False, scientific_eligibility=False, gpu=False,
        scope='Actual complete V6 software/job/source preparation; not native execution or completed '
              'startup admission. Uses frozen qualified future definitions; artificial copied job '
              'metadata is not misclassified as historical native seed usage. No original restart.')
    write_json(args.receipt, result)
    print(json.dumps({key:value for key,value in result.items() if key!='source_hashes'},indent=2))


if __name__ == '__main__':
    main()
