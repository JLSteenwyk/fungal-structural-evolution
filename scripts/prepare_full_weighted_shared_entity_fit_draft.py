#!/usr/bin/env python3
"""Prepare the complete unlaunched weighted fitting draft after software checks.

Real numerical closure, qualified-input timing and installed resources still
gate any launch. This script never launches a worker or modifies old plans.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from full_weighted_fit_exports import atomic
from full_weighted_covariance_sources_v2 import POLICIES, MODES
from reference_measurement_union_sources import verify


def run(output):
    assert not output.exists()
    gate = Path('metadata/full_weighted_fitting_software_validation_20261004_v2.json')
    proof = Path('metadata/full_weighted_fitting_software_transport_20261004_v2.json')
    validation = json.loads(gate.read_text()); transport = json.loads(proof.read_text())
    assert validation['status'] == 'passed_complete_four_control_fit_exports_readback_and_checkpoint_contracts_v2'
    assert validation['total_candidate_rows'] == 72000 and validation['total_setting_fit_links'] == 144000
    assert validation['actual_serialized_numeric_replays'] == 64
    assert validation['complete_export_grids_native_producer_and_reader_explicitly_mocked'] is True
    assert validation['production_fitting_launched'] is validation['scientific_eligibility'] is False
    assert transport['actual_tool_terminal_exit_code'] == 0
    verify(validation['source_hashes']); verify(transport['source_hashes'])
    source = Path('metadata/full_weighted_covariance_qualification_plan_20261004_v1.json')
    q = json.loads(source.read_text()); verify(q['pins'])
    old = Path('metadata/full_retained_shared_entity_fit_draft_plan_20261003_v1.json')
    production = json.loads(old.read_text()); assert production['launch_state'] == 'not_launched_or_queued'
    assert production['optimizer']['gradient_tolerance'] == 1e-6 and production['independent_audit']['column_batch'] == 32
    inventory = Path('metadata/full_weighted_shared_entity_fit_resource_inventory_20261004_v1.json')
    resources = json.loads(inventory.read_text())
    assert resources['candidate_rows'] == 20832000 and resources['setting_fit_links'] == 49766400
    assert resources['fits_launched_or_queued'] is resources['fit_resources_installed'] is False
    modules = ['full_weighted_fit_exports', 'full_weighted_shared_entity_fit_sources', 'prepare_full_weighted_shared_entity_fits',
        'readback_full_weighted_shared_entity_fits', 'weighted_shared_entity_candidate', 'readback_weighted_shared_entity_candidate',
        'shared_entity_likelihood', 'fit_shared_entity_likelihood', 'independent_shared_entity_likelihood',
        'independent_shared_entity_likelihood_fast', 'independent_shared_entity_optimizer', 'prepare_full_weighted_shared_entity_fit_draft',
        'reference_measurement_union_sources', 'full_weighted_covariance_qualification', 'full_weighted_covariance_sources_v2',
        'positive_diagonal_basis_context', 'independent_positive_diagonal_basis_context', 'readback_full_covariance_qualification']
    environment = Path('environments/full-shared-entity-fits-20261002.yml')
    paths = [gate, proof, source, old, inventory, environment, *[Path('scripts/' + m + '.py') for m in modules]]
    plan = dict(status='software_qualified_complete_weighted_fitting_draft_pending_real_closure_and_timing',
        prepared_utc=datetime.now(timezone.utc).isoformat(), qualification_plan=str(source),
        qualification_completion='metadata/full_weighted_covariance_qualification_completed_20261004_v1.json',
        expected=dict(logical_cases=75188, cohorts=4340, designs=130200, fit_inputs=260400, settings=622080,
            candidate_rows=20832000, setting_fit_links=49766400), methods=['ml','reml'], policies=POLICIES,
        loading_modes=MODES, trees=q['trees'], optimizer=production['optimizer'], independent_audit=production['independent_audit'],
        independent_backend='component_spectral_v1', output='results/phylogeny/full-four-control-shared-entity-fits-20261004-v1',
        environment=str(environment), fixture_validation=str(gate), source_production_configuration=str(old),
        resources=dict(cpus_per_worker=2,memory_gib_per_worker=32,swap_gib=0,blas_threads=1,address_space_gib=24,
            output_scratch_reserve_gib=3072,minimum_free_disk_gib=3172,runtime_calibrated=False,resource_inventory=str(inventory)),
        pins={str(p):sha(p) for p in paths}, launch_state='not_launched_or_queued', fits_computed=0,
        full_timing_implemented=False, scientific_eligibility=False, nonuniform_weighting_accepted=False,
        component_variance_attribution_accepted=False, gpu=False,new_cost_usd=0,
        scope='Complete original4340cohorts/130200designs/260400response inputs/622080settings, allfour policies/twomodes/fivetrees/twooutcomes andML/REML:20832000potential candidates and49766400original-setting links. All nonfits/reviews/failures retained, actual D/source audit/exact certificates bound, complete segmented exports and independent spectral/start/curvature/search/link replay software qualified. Source/journal fixtures and complete-grid native mocks are synthetic;64actual saved numerical fits separately replayed. Existing stricter production gradient1e-6/batch32 preserved. Real full numerical closure, complete qualified-input timing and installed resources remain mandatory before any native fitting launch. No worker queued, no fit resource installed, no ETA or accepted effect; all eight biological aims remain required.')
    assert not Path(plan['output']).exists()
    atomic(output, plan)
    print(json.dumps({k:v for k,v in plan.items() if k not in ['pins','scope']},indent=2))
    return plan


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True); run(p.parse_args().output)
