#!/usr/bin/env python3
"""Full V6 startup design, three actual native checks, explicit synthetic serialization."""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys
from unittest.mock import patch

from ancestral_chain_attempt import run_attempt, sha, write_json
from baliphy_scalar_json_logger_v6c import SCHEMA, restore
from check_baliphy_scalar_json_logger_v6_v9 import seed_census
from prepare_baliphy_scalar_v6_preflight_v1 import build, project_sources
from reference_measurement_union_sources import bind, verify
import run_baliphy_scalar_v6_preflight_v1 as producer
import readback_baliphy_scalar_v6_preflight_v1 as reader


def rejected(action):
    try: action()
    except (AssertionError, KeyError, ValueError, FileNotFoundError, TypeError): return
    raise AssertionError('Invalid startup state accepted')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args(); assert not a.receipt.exists()
    root = a.output.resolve(); root.mkdir(exist_ok=False)
    future_path = Path('metadata/baliphy_scalar_json_v6_future_models_20261004_v1.json')
    previous_path = Path('metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json')
    future = json.loads(future_path.read_text()); previous = json.loads(previous_path.read_text())
    pins = {}
    for mapping in [future['pins'], previous['pins']]:
        verify(mapping)
        for name, digest in mapping.items(): bind(pins, name, digest)
    modules = project_sources(pins, [Path(__file__),
        Path('scripts/prepare_baliphy_scalar_v6_preflight_v1.py'),
        Path('scripts/run_baliphy_scalar_v6_preflight_v1.py'),
        Path('scripts/readback_baliphy_scalar_v6_preflight_v1.py')])
    roles = json.loads(Path(future['future_roles']).read_text())
    originals = json.loads(Path(previous['jobs']).read_text())
    forbidden = seed_census(pins); assert len(forbidden) == future['forbidden_seed_count']
    jobs = build(roles, originals, forbidden)
    job_path = root / 'actual_full_grid_jobs.json'; write_json(job_path, jobs)
    altered = []
    for case in ['missing_role', 'duplicate_id', 'duplicate_source', 'duplicate_seed',
                 'old_seed', 'bool_seed', 'changed_prior', 'changed_alias', 'changed_role',
                 'wrong_seed_namespace', 'changed_program_digest', 'claimed_native_execution']:
        changed = copy.deepcopy(roles)
        if case == 'missing_role': changed.pop()
        elif case == 'duplicate_id': changed[1]['chain']['chain_id'] = changed[0]['chain']['chain_id']
        elif case == 'duplicate_source': changed[1]['source_v5_chain_id'] = changed[0]['source_v5_chain_id']
        elif case == 'duplicate_seed': changed[1]['chain']['seed'] = changed[0]['chain']['seed']
        elif case == 'old_seed': changed[0]['chain']['seed'] = originals[0]['chain']['seed']
        elif case == 'bool_seed': changed[0]['chain']['seed'] = True
        elif case == 'changed_prior': changed[0]['chain']['prior_label'] = 'invented'
        elif case == 'changed_alias': changed[0]['chain']['original_configuration_ids'] = []
        elif case == 'changed_role': changed[0]['chain']['chain'] = 5
        elif case == 'wrong_seed_namespace': changed[0]['seed_namespace'] = 'old'
        elif case == 'changed_program_digest': changed[0]['chain']['program_sha256'] = '0' * 64
        else: changed[0]['native_execution_launched'] = True
        rejected(lambda: build(changed, originals, forbidden)); altered.append(case)
    gate_path = Path('metadata/baliphy_scalar_json_v6_software_validation_20261004_v9.json')
    gate = json.loads(gate_path.read_text()); verify(gate['source_hashes'])
    for name, digest in gate['source_hashes'].items(): bind(pins, name, digest)
    fixtures = Path('data/software_audits/baliphy-joint-node-logger-20261003-v4').resolve()
    baseline = Path('data/software_audits/baliphy-joint-node-logger-20261003-v5').resolve()
    binary = Path(originals[0]['config']['command'][5])
    api = binary.parent.parent / 'lib/bali-phy/haskell'
    fixture_seeds = [20269011, 20269012, 20269013]
    assert set(fixture_seeds).isdisjoint(forbidden | {r['chain']['seed'] for r in roles})
    actual = {}; native_jobs = {}; native_receipts = {}; native_configs = {}
    for index, item in enumerate(gate['paired_prior_checks']):
        prior = item['prior']; directory = Path(item['new_directory'])
        prior_config = json.loads((directory.parent.parent / 'configuration.json').read_text())
        program = prior_config['command'][prior_config['command'].index('run') + 1]
        old = baseline / (prior + '-numeric-only.hs')
        assert restore(Path(program).read_text()) == old.read_text()
        paths = [Path('/usr/bin/prlimit'), binary, Path(program), old,
                 fixtures / 'alignment.faa', fixtures / 'tree.nwk', *sorted(api.rglob('*.hs'))]
        config = dict(command=['/usr/bin/prlimit', '--as=' + str(12 * 2**30),
            '--cpu=120', '--fsize=' + str(64 * 2**20), '--', str(binary),
            '--seed', str(fixture_seeds[index]), 'run', program, '--test', '--log-format', 'json'],
            timeout_seconds=180, pins={str(path): sha(path) for path in paths})
        chain = dict(chain_id='software-v6-startup-' + prior, effective_input_group='software',
            chain=1, prior_label=prior, family='synthetic', proteins=5,
            original_configuration_ids=['software-' + prior], seed=20267001 + index,
            program=str(old), program_sha256=sha(old), alignment=str(fixtures / 'alignment.faa'),
            alignment_sha256=sha(fixtures / 'alignment.faa'), tree=str(fixtures / 'tree.nwk'),
            tree_sha256=sha(fixtures / 'tree.nwk'))
        job = dict(chain=chain, fresh_seed=fixture_seeds[index], program=program,
                   project_scalar_schema=SCHEMA, config=config)
        receipt = run_attempt(root / 'native' / prior, config)
        row = producer.inspect(job, receipt, 'software-v6-startup')
        assert row['status'] == 'reference_startup_homology_density_and_representation_checked', row
        assert row['audit']['reference_homology_preserved'] and row['audit']['ancestral_lengths_free'] == 4
        actual[prior] = row; native_jobs[prior] = job
        native_receipts[prior] = receipt; native_configs[prior] = config
        for path in paths: bind(pins, path)
        print('actual_v6_native_startup', prior, row['status'], flush=True)
    output_cases = []
    for case in ['duplicate_key', 'bare_infinity', 'overflow_number', 'changed_density',
                 'unknown_footer', 'wrong_iteration', 'synthetic_native_failure']:
        folder = root / 'private_output_cases' / case
        original = native_receipts['broad']; shutil.copytree(original.parent.parent, folder)
        receipt = folder / 'attempt-0001' / 'receipt.json'
        stdout = receipt.parent / 'stdout.log'; text = stdout.read_text()
        value, end = json.JSONDecoder().raw_decode(text.lstrip()); footer = text.lstrip()[end:]
        if case == 'duplicate_key': text = json.dumps(value).replace('"iter": 0', '"iter": 0, "iter": 0', 1) + footer
        elif case in ['bare_infinity', 'overflow_number']:
            value['software_probe'] = 'replace-number'; text = json.dumps(value).replace('"replace-number"', 'Infinity' if case == 'bare_infinity' else '1e1000') + footer
        elif case == 'changed_density': value['parameters/']['referenceInitialization/']['priorFromDirectFormula'] += 1; text = json.dumps(value) + footer
        elif case == 'unknown_footer': text += '\nUnexpected footer\n'
        elif case == 'wrong_iteration': value['iter'] = 1; text = json.dumps(value) + footer
        stdout.write_text(text)
        native = json.loads(receipt.read_text())
        native['artifacts'] = {str(path.relative_to(receipt.parent)): sha(path)
            for path in receipt.parent.rglob('*') if path.is_file() and path != receipt}
        if case == 'synthetic_native_failure': native.update(exit_code=1, status='failed')
        write_json(receipt, native)
        row = producer.inspect(native_jobs['broad'], receipt, 'software-v6-startup')
        expected = 'unsuccessful_native_startup_retained' if case == 'synthetic_native_failure' else 'invalid_native_startup_retained'
        assert row['status'] == expected and row['scientific_eligibility'] is False, row
        output_cases.append(dict(case=case, status=row['status']))
    custody_cases = []
    for case in ['changed_configuration', 'changed_native_command', 'changed_native_process', 'unhashed_stdout']:
        folder = root / 'private_custody_cases' / case
        original = native_receipts['broad']; shutil.copytree(original.parent.parent, folder)
        receipt = folder / 'attempt-0001' / 'receipt.json'
        if case == 'unhashed_stdout': (receipt.parent / 'stdout.log').write_text('changed')
        else:
            path = folder / 'configuration.json' if case == 'changed_configuration' else receipt.parent / ('command.json' if case == 'changed_native_command' else 'process.json')
            value = json.loads(path.read_text())
            command = value if isinstance(value, list) else value['command']
            command[command.index('--seed') + 1] = '1'; write_json(path, value)
            native = json.loads(receipt.read_text())
            native['artifacts'] = {str(path.relative_to(receipt.parent)): sha(path)
                for path in receipt.parent.rglob('*') if path.is_file() and path != receipt}
            write_json(receipt, native)
        rejected(lambda: producer.inspect(native_jobs['broad'], receipt, 'software-v6-startup'))
        custody_cases.append(case)
    # Full metadata layout and actual V6 job construction are checked above.
    # The following controller/reader test explicitly mocks per-role native
    # execution and admission, and never represents full native startup.
    mocked_jobs = copy.deepcopy(jobs); fake = {}
    for index, job in enumerate(mocked_jobs):
        c = job['chain']; prior = c['prior_label']; job['config'] = native_configs[prior]
        row = copy.deepcopy(actual[prior]); row.update(chain_id=c['chain_id'],
            effective_input_group=c['effective_input_group'],
            model_input_identity=c['effective_input_group'] + '-' + prior,
            chain_role=c['chain'], family=c['family'], original_configuration_ids=c['original_configuration_ids'],
            source_seed=c['seed'], fresh_seed=job['fresh_seed'])
        if index in [0, 4]: row.update(status='unsuccessful_native_startup_retained', exit_code=1)
        fake[c['chain_id']] = row
    stage = root / 'explicitly_synthetic_full_stage'; synthetic_jobs = root / 'synthetic_jobs.json'
    write_json(synthetic_jobs, mocked_jobs)
    plan_path = root / 'synthetic_plan.json'
    write_json(plan_path, dict(jobs=str(synthetic_jobs), output=str(stage), pins={},
        resources=dict(memory_gib=16, minimum_free_disk_gib=1),
        scope='Explicitly mocked per-role native execution and source admission only.'))
    def fake_inspect(job, receipt, digest):
        row = copy.deepcopy(fake[job['chain']['chain_id']]); row['plan_sha256'] = digest; return row
    serial = []
    with patch.object(producer, 'run_attempt', side_effect=lambda root, config: native_receipts[next(k for k,v in native_configs.items() if v == config)]), \
         patch.object(producer, 'inspect', side_effect=fake_inspect), \
         patch.object(reader, 'inspect', side_effect=fake_inspect), \
         patch.object(sys, 'argv', ['software', '--plan', str(plan_path)]):
        producer.main(); summary = json.loads((stage / 'receipt.json').read_text())
        assert summary['validated_startups'] == 1618 and summary['unresolved_startup_quartets'] == 2
        rejected(producer.main); serial.append('producer_restart')
        for case in ['missing_checkpoint', 'changed_seed', 'scientific_acceptance', 'false_posterior']:
            cp = next((stage / 'chains').glob('*.json')); original = cp.read_bytes()
            if case == 'missing_checkpoint': cp.unlink()
            else:
                row = json.loads(original)
                if case == 'changed_seed': row['fresh_seed'] += 1
                elif case == 'scientific_acceptance': row['scientific_eligibility'] = True
                else: row['posterior_sampling_launched'] = True
                write_json(cp, row)
            rejected(lambda: reader.readback(plan_path)); cp.write_bytes(original); serial.append(case)
        reread = reader.readback(plan_path); assert reread['validated_startups'] == 1618
        rejected(lambda: reader.readback(plan_path)); serial.append('reader_restart')
    for path in [future_path, previous_path, Path(future['future_roles']), Path(previous['jobs']), gate_path, job_path]: bind(pins, path)
    for path in root.rglob('*'):
        if path.is_file(): bind(pins, path)
    verify(pins)
    result = dict(status='passed_full_scalar_v6_startup_software_contracts_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(), full_roles=1620, full_quartets=405,
        effective_inputs=135, original_configuration_aliases=324,
        full_source_grid_checked=True, native_startups=3, native_fixture_seeds=fixture_seeds,
        actual_native_startup_checks=list(actual.values()), forbidden_seeds=sorted(forbidden),
        forbidden_seed_count=len(forbidden), transitive_project_source_modules=len(modules),
        altered_designs_rejected=altered, privately_rehashed_invalid_output_cases=output_cases,
        custody_cases_rejected=custody_cases, synthetic_full_producer_reader_serialization_checked=True,
        artificial_failures_retained=2, artificial_unresolved_quartets_retained=2,
        serialization_rejections=serial, source_hashes=pins,
        full_grid_native_startup_launched=False, scientific_eligibility=False,
        posterior_sampling_launched=False, gpu=False, new_cost_usd=0,
        scope='All1620real prepared V6 startup role constructions. Three actual native--test '
              'synthetic fixtures across priors validate preserved initialization/densities/representation '
              'and rooted trees. Seven rehashed private invalid/failure outputs retained and four '
              'custody corruptions rejected. Full1620controller/reader serialization explicitly mocks '
              'native execution and per-role admission, retaining2failures/2unresolvedquartets. '
              'Not a fungal pilot, fullnativegrid, MCMC scalar/joint logging check, posterior '
              'qualification, original restart or biological acceptance.')
    with a.receipt.open('x') as f: json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes', 'forbidden_seeds', 'actual_native_startup_checks']}, indent=2))


if __name__ == '__main__':
    main()
