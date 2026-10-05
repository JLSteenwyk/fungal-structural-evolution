#!/usr/bin/env python3
"""Read all current short-run alpha traces and execute a deterministic native probe."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys

from ancestral_chain_attempt import run_attempt, sha
from read_baliphy_scalar_json_v6b import load, read_record
from reference_measurement_union_sources import bind, verify


PRIORS = {'package': (6., 2.), 'centered': (0., 1.), 'broad': (0., 2.)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists() and not args.receipt.exists()
    args.output.mkdir(parents=True)
    pins = {}
    source_specs = [
        ('v6', Path('metadata/baliphy_scalar_v6_sampler_execution_plan_20261004_v1.json'),
         Path('metadata/baliphy_scalar_v6_short_sampler_completed_20261004_v1.json'), 1620),
        ('v7_comparison', Path('metadata/baliphy_joint_fasta_v7_failure_grid_plan_20261004_v1.json'),
         Path('metadata/baliphy_joint_fasta_v7_failure_grid_completed_20261005_v1.json'), 24)]
    counts = Counter()
    per_role = []
    special = []
    expected_header = {'fields': ['iter','prior','likelihood','posterior'], 'nested': True,
        'format':'MCON', 'version':'0.2', 'projectScalarSchema':'native-cjson-explicit-special-values-v6'}
    trace = args.output/'all_current_alpha_traces.jsonl'
    with trace.open('x') as handle:
        for stage, plan_path, completion_path, expected_roles in source_specs:
            plan = json.loads(plan_path.read_text())
            completion = json.loads(completion_path.read_text())
            assert completion['scientific_eligibility'] is False
            if stage == 'v6':
                assert completion['status'] == 'complete_verified_full_scalar_v6_short_sampler_qualification_v1'
                assert sha(completion['full_hash_archive']) == completion['full_hash_archive_sha256']
                bind(pins, completion['full_hash_archive'], completion['full_hash_archive_sha256'])
            else:
                assert completion['status'] == 'complete_verified_all24_v7_failure_comparison_dispositions'
            jobs_path = Path(plan['jobs'])
            jobs = json.loads(jobs_path.read_text())
            jobs = {item['chain']['chain_id']: item for item in jobs}
            disposition_path = Path(plan['output'])/'dispositions.json'
            dispositions = json.loads(disposition_path.read_text())
            assert len(dispositions) == len(jobs) == expected_roles
            assert {row['chain_id'] for row in dispositions} == set(jobs)
            for path in [plan_path, completion_path, jobs_path, disposition_path]:
                bind(pins, path)
            for role in dispositions:
                job = jobs[role['chain_id']]
                chain = job['chain']
                program = Path(chain['program'])
                bind(pins, program, chain['program_sha256'])
                mu, scale = PRIORS[role['prior_label']]
                # The actual generated prior is the continuous exp transform,
                # not a spike at infinity or a changed/truncated prior.
                source = program.read_text()
                expected = ';alpha_2 <- sample (logLaplace '+str(int(mu))+' '+str(int(scale))+')'
                assert source.count(expected) == 1
                receipt_path = Path(role['native_receipt'])
                bind(pins, receipt_path, role['native_receipt_sha256'])
                native = json.loads(receipt_path.read_text())
                assert native['exit_code'] == role['exit_code']
                counts[stage+'_roles'] += 1
                counts['native_zero_roles' if native['exit_code']==0 else 'retained_native_failures'] += 1
                row_summary = dict(stage=stage, chain_id=role['chain_id'], prior=role['prior_label'],
                    native_exit_code=native['exit_code'], original_disposition=role['status'],
                    scalar_integrity_accepted=role['scalar_integrity_accepted'], rows=0,
                    nonfinite_alpha_iterations=[], maximum_finite_alpha=None)
                directories = list(receipt_path.parent.glob('independent-chain-*'))
                assert len(directories) == 1
                directory = directories[0]
                jp, tp, mp = [directory/name for name in ['C1.log.json','C1.log','C1.log.column-map.json']]
                for path in [jp, tp, mp]:
                    if path.exists():
                        relative = str(path.relative_to(receipt_path.parent))
                        bind(pins, path, native['artifacts'][relative])
                if not all(path.exists() for path in [jp,tp,mp]):
                    assert native['exit_code'] != 0
                    row_summary['unavailable_scalar_files'] = [str(path) for path in [jp,tp,mp] if not path.exists()]
                    per_role.append(row_summary)
                    continue
                lines = jp.read_text().splitlines()
                assert load(lines[0]) == expected_header
                records = [load(line) for line in lines[1:]]
                with tp.open() as tsv_handle:
                    tsv = list(csv.DictReader(tsv_handle, delimiter='\t'))
                mapping = load(mp.read_text())
                assert len(records) == len(tsv)
                iterations = [record['iter'] for record in records]
                assert iterations == list(range(len(records)))
                if native['exit_code'] == 0:
                    assert iterations == list(range(21))
                finite = []
                for record, reference in zip(records, tsv):
                    audit = read_record(record)
                    iteration = record['iter']
                    assert int(reference[mapping['iter']]) == iteration
                    alpha = record['parameters//']['S1/']['ASRV.Gamma:alpha']
                    native_alpha_token = reference[mapping['S1/ASRV.Gamma:alpha']]
                    native_alpha = float(native_alpha_token)
                    reviews = audit['nonfinite_reviews']
                    assert all(item == dict(section='parameters', path=['S1/','ASRV.Gamma:alpha'],
                                            kind='positive_infinity') for item in reviews)
                    if reviews:
                        assert len(reviews)==1 and alpha=='__project_scalar_v6__:positive_infinity'
                        assert math.isinf(native_alpha) and native_alpha > 0
                        alpha_value = None
                        counts['positive_infinity_alpha_observations'] += 1
                        row_summary['nonfinite_alpha_iterations'].append(iteration)
                    else:
                        assert type(alpha) in (int,float) and math.isfinite(alpha) and alpha > 0
                        assert math.isclose(alpha,native_alpha,rel_tol=2e-13,abs_tol=0)
                        alpha_value = alpha
                        finite.append(alpha)
                    context = record['statistics//']
                    logs = {}
                    for key in ['prior','likelihood','posterior']:
                        value = context[key]
                        assert type(value) in (int,float) and math.isfinite(value)
                        expected_value = float(reference[mapping[key]])
                        assert math.isfinite(expected_value) and math.isclose(value,expected_value,rel_tol=2e-13,abs_tol=0)
                        logs[key] = value
                    assert math.isclose(logs['posterior'],logs['prior']+logs['likelihood'],rel_tol=2e-13,abs_tol=1e-8)
                    emitted = dict(stage=stage,chain_id=role['chain_id'],family=role['family'],
                        prior_label=role['prior_label'],iteration=iteration,
                        alpha=alpha_value,alpha_state='positive_infinity' if reviews else 'finite',
                        native_alpha_token=native_alpha_token,
                        reconstructed_log_alpha=math.log(alpha_value) if alpha_value is not None else None,
                        latent_log_alpha_observed=False,**logs)
                    handle.write(json.dumps(emitted,allow_nan=False)+'\n')
                    counts['scalar_rows'] += 1
                    if reviews:special.append(emitted)
                row_summary['rows'] = len(records)
                row_summary['maximum_finite_alpha'] = max(finite) if finite else None
                per_role.append(row_summary)
    assert counts['v6_roles']==1620 and counts['v7_comparison_roles']==24
    assert counts['native_zero_roles']==1620 and counts['retained_native_failures']==24
    assert counts['scalar_rows']==34020 and counts['positive_infinity_alpha_observations']==26
    assert sum(bool(role['nonfinite_alpha_iterations']) for role in per_role)==14
    binary = Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/bali-phy').resolve()
    template = Path('config/native-format-probes/log-alpha-overflow-v1.hs')
    program = (args.output/'log-alpha-overflow-v1.hs').resolve()
    program.write_bytes(template.read_bytes())
    library = binary.parent.parent/'lib/bali-phy'
    native_paths = [Path('/usr/bin/prlimit'),binary,program,*library.glob('*.so'),
        *[library/'haskell'/path for path in ['SModel/ASRV.hs','Probability/Distribution/Laplace.hs',
             'Probability/Distribution/Transform.hs','Numeric/Log.hs']]]
    config = dict(command=['/usr/bin/prlimit','--as='+str(8*2**30),'--cpu=60','--fsize='+str(16*2**20),
        '--',str(binary),'run',str(program)],timeout_seconds=90,pins={str(path):sha(path) for path in native_paths})
    native_receipt = run_attempt(args.output/'deterministic-native-probe',config)
    native = json.loads(native_receipt.read_text())
    assert native['exit_code']==0
    output = native_receipt.parent/'stdout.log'
    observed = [load(line) for line in output.read_text().splitlines()]
    assert len(observed)==57
    grid = [-20.,-10.,-6.,-3.,0.,6.,20.,50.,100.,500.,700.,709.,709.7,709.78,709.79,710.,750.,1000.,5000.]
    native_infinite = 0
    for label,(mu,scale) in PRIORS.items():
        subset = [row for row in observed if row['prior']==label]
        assert [row['log_alpha'] for row in subset]==grid
        for row in subset:
            x = row['log_alpha']
            expected_density = -math.log(2*scale)-abs(x-mu)/scale
            assert row['latent_log_density_tag']=='finite'
            assert math.isclose(row['latent_log_density'],expected_density,rel_tol=2e-13,abs_tol=2e-13)
            try:expected_alpha=math.exp(x)
            except OverflowError:expected_alpha=math.inf
            if math.isinf(expected_alpha):
                assert row['alpha_tag']=='positive_infinity' and row['alpha'] is None
                assert row['derived_log_density_tag']=='negative_infinity' and row['derived_log_density'] is None
                native_infinite += 1
            else:
                assert row['alpha_tag']=='finite'
                assert math.isclose(row['alpha'],expected_alpha,rel_tol=2e-13,abs_tol=0)
                expected_pushforward = expected_density-x
                assert row['derived_log_density_tag']=='finite'
                assert math.isclose(row['derived_log_density'],expected_pushforward,rel_tol=2e-13,abs_tol=2e-13)
            rates = row['gamma_rates']
            assert len(rates)==4 and all(type(v) in (int,float) and math.isfinite(v) and v>=0 for v in rates)
            assert math.isclose(math.fsum(rates)/4,1.,rel_tol=1e-10,abs_tol=1e-10)
            if x>=100:assert rates==[1.,1.,1.,1.]
    assert native_infinite==15
    for path,digest in config['pins'].items():bind(pins,path,digest)
    bind(pins,native_receipt)
    for path,digest in native['artifacts'].items():bind(pins,native_receipt.parent/path,digest)
    for path in [Path(__file__),Path('scripts/read_baliphy_scalar_json_v6b.py'),template,trace]:bind(pins,path)
    verify(pins)
    result = dict(status='completed_full_current_alpha_trace_and_native_overflow_diagnostic',
        checked_utc=datetime.now(timezone.utc).isoformat(),counts=dict(counts),role_summaries=per_role,
        special_observations=special,native_deterministic_cases=57,native_exp_overflow_cases=15,
        continuous_generated_priors_checked=411,source_hashes=pins,scientific_eligibility=False,
        latent_log_alpha_recorded_in_original_samples=False,original_chain_overflow_cause_proven=False,
        model_or_prior_changed=False,reviewed_arrays_admitted=False,new_mcmc_runs=0,gpu=False,
        scope='All1620V6roles and24V7comparison roles retained, including24original native failures. '
              'Every available scalar trace and actual generated prior bound; JSON/TSV alpha and finite '
              'prior/likelihood/posterior checks preserve original tolerances.57deterministic native cases '
              'show finite latent Laplace density can coexist with exp overflow and unit GammaMean rates. '
              'Original latent log-alpha was not saved: this identifies a compatible mechanism, not its '
              'historical causal proof, a point mass at infinity, an installed formatter repair, posterior '
              'adequacy, scalar review admission or completed biological aim.')
    with args.receipt.open('x') as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','role_summaries','special_observations']},indent=2))


if __name__=='__main__':
    main()
