#!/usr/bin/env python3
"""Qualify complete V6 sampler serialization with retained actual fixtures and mocks."""
import argparse
import copy
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

from ancestral_chain_attempt import sha, write_json
from baliphy_joint_sampler_scalar_v6 import build_jobs, inspect, SUCCESS, INVALID, REVIEW
from baliphy_scalar_json_logger_v6c import SCHEMA
from baliphy_scalar_v6_sampler_gates_v1 import prerequisites, startup_roles
from check_baliphy_scalar_json_logger_v6_v9 import seed_census
from independent_joint_ancestral_frames import verify_arrays, write_arrays
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind, verify
from reference_sampler_memory_budget import MemoryBudget
import run_baliphy_scalar_v6_sampler_v1 as workflow


def rejected(action):
    try: action()
    except (AssertionError, KeyError, ValueError, RuntimeError, FileNotFoundError): return
    raise AssertionError('Altered controller state accepted')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args(); assert not a.receipt.exists()
    root = a.output.resolve(); root.mkdir(exist_ok=False)
    future_path = Path('metadata/baliphy_scalar_json_v6_future_models_20261004_v1.json')
    previous_path = Path('metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json')
    admission_path = Path('metadata/baliphy_scalar_v6_joint_adapter_software_transport_20261004_v2.json')
    future = json.loads(future_path.read_text()); previous = json.loads(previous_path.read_text())
    admission = json.loads(admission_path.read_text()); pins = {}
    assert admission['status'] == 'verified_original_v6_scalar_joint_adapter_v2_software_wait_zero'
    for mapping in [future['pins'], previous['pins'], admission['source_hashes']]:
        verify(mapping)
        for name, digest in mapping.items(): bind(pins, name, digest)
    modules = project_sources(pins, [Path(__file__), Path('scripts/run_baliphy_scalar_v6_sampler_v1.py'),
        Path('scripts/baliphy_scalar_v6_sampler_gates_v1.py')])
    roles = json.loads(Path(future['future_roles']).read_text())
    originals = json.loads(Path(previous['jobs']).read_text())
    forbidden = seed_census(pins); assert len(forbidden) == future['forbidden_seed_count']
    jobs = build_jobs(roles, originals, forbidden)
    job_path = root / 'actual_full_grid_jobs.json'; write_json(job_path, jobs)
    original_checked_jobs = Path('data/software_audits/baliphy-scalar-json-v6-joint-adapter-20261004-v2/software_grid_jobs.json').resolve()
    assert sha(job_path) == sha(original_checked_jobs)
    mapping = root / 'fixture_mapping.tsv'
    with mapping.open('x') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=['guide','family','dataset','level','source_node','retained_set_json'])
        writer.writeheader()
        for i,n in enumerate([2,3,4,5]): writer.writerow(dict(guide='profile', family='software', dataset='whole',
            level=i, source_node='n'+str(i), retained_set_json=json.dumps(list('abcde')[:n])))
    fixture = Path('data/software_audits/baliphy-joint-node-logger-20261003-v4').resolve()
    scalar_gate_path = Path('metadata/baliphy_scalar_json_v6_software_validation_20261004_v9.json')
    scalar_gate = json.loads(scalar_gate_path.read_text())
    actual = {}; actual_jobs = {}; receipts = {}; arrays = {}
    for item in scalar_gate['paired_prior_checks']:
        prior = item['prior']; directory = Path(item['new_directory']); receipt = directory.parent / 'receipt.json'
        config = json.loads((receipt.parent.parent / 'configuration.json').read_text())
        program = config['command'][config['command'].index('run')+1]
        base = next(j['chain'] for j in jobs if j['chain']['prior_label'] == prior)
        chain = dict(base, chain_id='software-v6-controller-'+prior, seed=item['seed'],
            family='software', proteins=5, effective_input_group='software',
            original_configuration_ids=['software-whole-'+prior], program=program, program_sha256=sha(program),
            alignment=str(fixture/'alignment.faa'), alignment_sha256=sha(fixture/'alignment.faa'),
            tree=str(fixture/'tree.nwk'), tree_sha256=sha(fixture/'tree.nwk'))
        job = dict(chain=chain, config=config, source_seed=20267001, memory_reservation_bytes=12*2**30)
        export = root / 'actual_fixture_arrays' / prior
        row = inspect(job, receipt, 'software-v6-controller', mapping, export)
        assert row['status'] == SUCCESS and row['scalar_integrity_accepted']
        assert row['scalar_v6_audit']['rows'] == 21 and row['scalar_v6_audit']['mapped_values_compared'] == 903
        assert inspect(job, receipt, 'software-v6-controller', mapping, export, False) == row
        actual[prior] = row; actual_jobs[prior] = job; receipts[prior] = receipt; arrays[prior] = []
        for frame in row['joint_frames']:
            with np.load(frame['projection_array'], allow_pickle=False) as data:
                arrays[prior].append({k:data[k].copy() for k in data.files})
        absent = root / 'missing_actual_exports' / prior
        rejected(lambda: inspect(job, receipt, 'software-v6-controller', mapping, absent, False))
        assert not absent.exists()
    negative_source_root = Path('data/software_audits/baliphy-scalar-json-v6-joint-adapter-20261004-v2').resolve()
    retained_negative_rows = {}
    for case, expected in [('changed_finite_value', INVALID), ('tagged_infinity_review', REVIEW),
                           ('synthetic_native_failure', 'unsuccessful_sampler_qualification_attempt_retained')]:
        receipt = negative_source_root / 'altered' / case / 'attempt-0001' / 'receipt.json'
        assert str(receipt) in admission['source_hashes']
        assert sha(receipt) == admission['source_hashes'][str(receipt)]
        export = root / 'retained_negative_exports' / case
        row = inspect(actual_jobs['broad'], receipt, 'software-v6-controller', mapping, export)
        assert row['status'] == expected and not row['scalar_integrity_accepted'] and not row['joint_frames']
        assert not export.exists(); retained_negative_rows[case] = row
    # Full real role metadata and jobs above; only per-role native execution,
    # admission, prerequisite admission and runtime caps are mocked below.
    # Successful arrays are genuinely written and reread in this new audit root.
    output = root / 'explicitly_synthetic_full_controller'; plan_path = root / 'synthetic_plan.json'
    write_json(plan_path, dict(jobs=str(job_path), output=str(output), mapping=str(mapping), pins={},
        resources=dict(cpus=16, memory_gib=200, reservation_capacity_gib=192,
            workers=16, minimum_free_disk_gib=1), scope='Explicit mocked native/admission/gates/caps; real serialized arrays.'))
    negative_ids = {jobs[0]['chain']['chain_id']:'changed_finite_value',
                    jobs[4]['chain']['chain_id']:'tagged_infinity_review',
                    jobs[8]['chain']['chain_id']:'synthetic_native_failure'}
    by_id = {j['chain']['chain_id']:j for j in jobs}
    def mock_attempt(folder, config): return receipts[by_id[Path(folder).name]['chain']['prior_label']]
    def mock_inspect(job, receipt, digest, mapping, export, allow_export_creation=True):
        chain = job['chain']; prior = chain['prior_label']; cid = chain['chain_id']
        template = retained_negative_rows[negative_ids[cid]] if cid in negative_ids else actual[prior]
        row = copy.deepcopy(template)
        row.update(chain_id=cid, effective_input_group=chain['effective_input_group'],
            model_input_identity=chain['effective_input_group']+'-'+prior, chain_role=chain['chain'],
            family=chain['family'], prior_label=prior, original_configuration_ids=chain['original_configuration_ids'],
            seed=chain['seed'], source_seed=job['source_seed'], memory_reservation_bytes=job['memory_reservation_bytes'],
            plan_sha256=digest)
        if cid in negative_ids:
            assert not Path(export).exists(); return row
        export = Path(export); exists = export.exists()
        assert exists or allow_export_creation
        if not exists: export.mkdir(parents=True)
        for index, frame in enumerate(row['joint_frames']):
            path = export / ('frame-'+str(frame['iteration'])+'.npz')
            if not exists: write_arrays(path, arrays[prior][index])
            verify_arrays(path, arrays[prior][index])
            frame.update(projection_array=str(path), projection_array_sha256=sha(path))
        return row
    serial = []
    with patch.object(workflow, 'startup_gate', return_value={}), \
         patch.object(workflow, 'runtime_caps', return_value={}), \
         patch.object(workflow, 'run_attempt', side_effect=mock_attempt), \
         patch.object(workflow, 'inspect', side_effect=mock_inspect):
        producer = workflow.run(plan_path); reader = workflow.run(plan_path, reader=True)
        assert producer['checked_sampler_attempts'] == 1617 and producer['unsuccessful_sampler_attempts'] == 3
        assert producer['complete_quartets'] == 402 and producer['unresolved_quartets'] == 3
        assert producer['scalar_review_roles'] == 1 and producer['scalar_v6_finite_checked_roles'] == 1617
        assert producer['scalar_v6_mapped_values_checked'] == 1617*903 and producer['joint_saved_frames'] == 4851
        rejected(lambda: workflow.run(plan_path)); serial.append('producer_restart')
        rejected(lambda: workflow.run(plan_path, reader=True)); serial.append('reader_restart')
        # Preserve successful exported bytes. Alterations below are confined
        # to explicitly synthetic artifacts, never parent native fixtures.
        reader_path = output / 'readback.json'; reader_bytes = reader_path.read_bytes(); reader_path.unlink()
        cp = output / 'chains' / (jobs[1]['chain']['chain_id']+'.json'); cp_bytes = cp.read_bytes()
        row = json.loads(cp_bytes); frame_path = Path(row['joint_frames'][0]['projection_array'])
        frame_bytes = frame_path.read_bytes()
        for case in ['missing_checkpoint', 'missing_frame', 'missing_export_directory', 'changed_seed',
                     'scalar_acceptance_claim', 'scalar_comparison_count', 'scientific_claim', 'foreign_export', 'changed_ledger']:
            if case == 'missing_checkpoint': cp.unlink()
            elif case == 'missing_frame': frame_path.unlink()
            elif case == 'missing_export_directory':
                preserved = frame_path.parent.with_name(frame_path.parent.name+'-preserved'); frame_path.parent.rename(preserved)
            elif case == 'foreign_export':
                foreign = frame_path.parent / 'foreign.npz'; foreign.write_bytes(frame_bytes)
            elif case == 'changed_ledger':
                ledger = output / 'memory_reservations.jsonl'; ledger_bytes = ledger.read_bytes()
                events = ledger.read_text().splitlines(); event = json.loads(events[0]); event['reserved_total'] += 1
                events[0] = json.dumps(event); ledger.write_text('\n'.join(events)+'\n')
            else:
                changed = json.loads(cp_bytes)
                if case == 'changed_seed': changed['seed'] += 1
                elif case == 'scalar_acceptance_claim': changed['scalar_integrity_accepted'] = False
                elif case == 'scalar_comparison_count': changed['scalar_v6_audit']['mapped_values_compared'] += 1
                else: changed['scientific_eligibility'] = True
                write_json(cp, changed)
            rejected(lambda: workflow.run(plan_path, reader=True))
            if case == 'missing_checkpoint': cp.write_bytes(cp_bytes)
            elif case == 'missing_frame': assert not frame_path.exists(); frame_path.write_bytes(frame_bytes)
            elif case == 'missing_export_directory': assert not frame_path.parent.exists(); preserved.rename(frame_path.parent)
            elif case == 'foreign_export': foreign.unlink()
            elif case == 'changed_ledger': ledger.write_bytes(ledger_bytes)
            else: cp.write_bytes(cp_bytes)
            serial.append(case)
        reader_path.write_bytes(reader_bytes)
    # Check new startup metadata admission using explicitly synthetic rows,
    # never fabricated native completion journals.
    startup_root = root / 'synthetic_startup'; startup_root.mkdir()
    startup_path = startup_root / 'dispositions.json'
    rows = [dict(chain_id=j['chain']['chain_id'], source_seed=j['source_seed'], fresh_seed=j['chain']['seed'],
        chain_role=j['chain']['chain'], prior_label=j['chain']['prior_label'],
        effective_input_group=j['chain']['effective_input_group'], original_configuration_ids=j['chain']['original_configuration_ids'],
        project_scalar_schema=SCHEMA, scientific_eligibility=False, posterior_sampling_launched=False,
        status='reference_startup_homology_density_and_representation_checked') for j in jobs]
    closed = dict(full_chains=1620, full_quartets=405, effective_inputs=135, original_configuration_aliases=324,
        validated_startups=1620, unsuccessful_startups=0, complete_startup_quartets=405,
        unresolved_startup_quartets=0, posterior_sampling_launched=False)
    write_json(startup_path, rows)
    startup_roles(dict(output=str(startup_root)), closed, {str(startup_path):sha(startup_path)}, jobs)
    gate_cases = []
    for case in ['missing_role','duplicate_role','wrong_fresh_seed','wrong_source_seed','wrong_role','wrong_prior',
                 'wrong_alias','wrong_schema','native_failure','scientific_claim','posterior_claim']:
        changed = copy.deepcopy(rows)
        if case == 'missing_role': changed.pop()
        elif case == 'duplicate_role': changed[1] = copy.deepcopy(changed[0])
        elif case == 'wrong_fresh_seed': changed[0]['fresh_seed'] += 1
        elif case == 'wrong_source_seed': changed[0]['source_seed'] += 1
        elif case == 'wrong_role': changed[0]['chain_role'] = 5
        elif case == 'wrong_prior': changed[0]['prior_label'] = 'invented'
        elif case == 'wrong_alias': changed[0]['original_configuration_ids'] = []
        elif case == 'wrong_schema': changed[0]['project_scalar_schema'] = 'old'
        elif case == 'native_failure': changed[0]['status'] = 'unsuccessful_native_startup_retained'
        elif case == 'scientific_claim': changed[0]['scientific_eligibility'] = True
        else: changed[0]['posterior_sampling_launched'] = True
        write_json(startup_path, changed)
        rejected(lambda: startup_roles(dict(output=str(startup_root)), closed,
            {str(startup_path):sha(startup_path)}, jobs)); gate_cases.append(case)
    write_json(startup_path, rows)
    actual_gate_plan = dict(startup_plan='metadata/baliphy_scalar_v6_preflight_plan_20261004_v1.json',
        jobs=str(job_path), previous_joint_sampler_plan=str(previous_path))
    try:
        actual_bindings = prerequisites(actual_gate_plan)
        for name, digest in actual_bindings.items(): bind(pins, name, digest)
        actual_gate = 'all_original_prerequisite_closures_verified'
    except FileNotFoundError as error:
        assert error.filename == 'metadata/baliphy_scalar_v6_preflight_completed_20261004_v1.json'
        actual_gate = 'actual_unclosed_v6_startup_refused_before_native_admission'
    # An exceptional or incomplete native attempt aborts memory admission;
    # no following role is allowed to call the native runner.
    abort_cases = []
    for case in ['native_runner_exception', 'incomplete_existing_attempt']:
        stage = root / 'private_abort_cases' / case; stage.mkdir(parents=True)
        job = jobs[0]; budget = MemoryBudget(192*2**30)
        if case == 'incomplete_existing_attempt':
            attempt = stage / 'attempts' / job['chain']['chain_id'] / 'attempt-0001'
            attempt.mkdir(parents=True); write_json(attempt/'process.json', dict(synthetic=True))
        with patch.object(workflow, 'run_attempt', side_effect=RuntimeError('Synthetic native runner exception')) as native:
            rejected(lambda: workflow.execute_job(job, stage, 'synthetic-abort', mapping, budget, 1))
            assert budget.aborted is not None and budget.used == 0
            count = native.call_count
            rejected(lambda: workflow.execute_job(jobs[1], stage, 'synthetic-abort', mapping, budget, 1))
            assert native.call_count == count == (1 if case == 'native_runner_exception' else 0)
        write_json(stage/'memory_events.json', budget.events); abort_cases.append(case)
    for path in [future_path, previous_path, admission_path, scalar_gate_path, Path(future['future_roles']),
                 Path(previous['jobs']), original_checked_jobs, job_path, mapping, Path(actual_gate_plan['startup_plan'])]: bind(pins, path)
    for path in root.rglob('*'):
        if path.is_file(): bind(pins, path)
    verify(pins)
    result = dict(status='passed_full_scalar_v6_sampler_controller_software_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(), full_roles=1620, full_quartets=405,
        effective_inputs=135, original_configuration_aliases=324,
        full_job_file_matches_qualified_admission_adapter=True, inherited_design_rejections=12,
        retained_actual_native_roles=3, retained_actual_scalar_rows=63, retained_actual_mapped_values=2709,
        retained_actual_joint_frames=9, missing_actual_export_readbacks_rejected=3,
        retained_synthetic_negative_templates=3, synthetic_full_controller_serialization_checked=True,
        synthetic_checked_roles=1617, synthetic_unresolved_roles=3, synthetic_complete_quartets=402,
        synthetic_unresolved_quartets=3, synthetic_scalar_review_roles=1, synthetic_joint_frames=4851,
        synthetic_reservation_ledger=producer['reservation_audit'], serialization_cases_rejected=serial,
        synthetic_startup_metadata_cases_rejected=gate_cases, actual_prerequisite_gate=actual_gate,
        abort_admission_cases_checked=abort_cases, transitive_project_source_modules=len(modules), source_hashes=pins,
        new_native_runs=0, full_grid_sampler_launched=False, native_resource_observer_qualified=False,
        posterior_qualified=False, scientific_eligibility=False, gpu=False, new_cost_usd=0,
        scope='Complete real1620job construction matches previous qualified V6 admission bytes. '
              'Three retained actual native fixtures reread63scalar rows/2709mapped values/9jointframes. '
              'Full producer/reader serialization explicitly mocks native execution/admission/gates/caps; '
              'actual synthetic array writes/readback and reservation ledger retain invalid/review/failure '
              'roles. No invented original journals, future native sampler launch, measured future '
              'resource proof, posterior acceptance, original restart or biological inference.')
    with a.receipt.open('x') as f: json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','synthetic_reservation_ledger']}, indent=2))


if __name__ == '__main__':
    main()
