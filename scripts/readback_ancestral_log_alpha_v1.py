#!/usr/bin/env python3
"""Independently reconstruct all alpha traces and check native probes with Decimal."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
import re


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def strict(text):
    def pairs(items):
        output={}
        for key,value in items:
            if key in output:raise ValueError('Duplicate key')
            output[key]=value
        return output
    def invalid(value):raise ValueError('Invalid JSON number '+value)
    return json.loads(text,object_pairs_hook=pairs,parse_constant=invalid)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args()
    assert not args.receipt.exists()
    pp=Path('metadata/ancestral_log_alpha_diagnostic_20261005_v2.json')
    tp=Path('metadata/ancestral_log_alpha_diagnostic_transport_20261005_v2.json')
    producer,transport=[strict(path.read_text()) for path in [pp,tp]]
    assert producer['status']=='completed_full_current_alpha_trace_and_native_overflow_diagnostic_v2'
    assert transport['validation_sha256']==sha(pp) and transport['original_tool_terminal_exit_code']==0
    pins=dict(producer['source_hashes'])
    for path,digest in pins.items():assert sha(path)==digest,path
    trace=Path('results/ancestral/full-current-log-alpha-diagnostic-20261005-v2/all_current_alpha_traces.jsonl')
    priors={'package':(6,2),'centered':(0,1),'broad':(0,2)}
    counts=Counter();special=[];roles=[];programs=set();maximum_log_difference=0.
    plans=[('v6','metadata/baliphy_scalar_v6_sampler_execution_plan_20261004_v1.json',1620),
           ('v7_comparison','metadata/baliphy_joint_fasta_v7_failure_grid_plan_20261004_v1.json',24)]
    with localcontext() as decimal_context,trace.open() as trace_handle:
        decimal_context.prec=90
        emitted=(strict(line) for line in trace_handle)
        for stage,plan_path,expected_roles in plans:
            plan=strict(Path(plan_path).read_text())
            source=strict((Path(plan['output'])/'dispositions.json').read_text())
            jobs={item['chain']['chain_id']:item['chain'] for item in strict(Path(plan['jobs']).read_text())}
            assert len(source)==len(jobs)==expected_roles
            for role in source:
                counts[stage+'_roles']+=1
                counts['native_zero_roles' if role['exit_code']==0 else 'retained_native_failures']+=1
                chain=jobs[role['chain_id']]
                program=Path(chain['program']);programs.add(str(program))
                assert sha(program)==chain['program_sha256']
                definitions=re.findall(r';alpha_2 <- sample \(logLaplace ([0-9.]+) ([0-9.]+)\)',program.read_text())
                assert len(definitions)==1 and tuple(map(float,definitions[0]))==priors[role['prior_label']]
                native_path=Path(role['native_receipt']);assert sha(native_path)==role['native_receipt_sha256']
                native=strict(native_path.read_text());assert native['exit_code']==role['exit_code']
                paths=list(native_path.parent.glob('independent-chain-*'));assert len(paths)==1
                directory=paths[0]
                jp,tp_raw,mp=[directory/name for name in ['C1.log.json','C1.log','C1.log.column-map.json']]
                summary=dict(stage=stage,chain_id=role['chain_id'],prior=role['prior_label'],
                    native_exit_code=native['exit_code'],original_disposition=role['status'],
                    scalar_integrity_accepted=role['scalar_integrity_accepted'],rows=0,
                    nonfinite_alpha_iterations=[],maximum_finite_alpha=None)
                if not all(path.exists() for path in [jp,tp_raw,mp]):
                    assert native['exit_code']!=0
                    summary['unavailable_scalar_files']=[str(path) for path in [jp,tp_raw,mp] if not path.exists()]
                    roles.append(summary);continue
                lines=jp.read_text().splitlines();header=strict(lines[0])
                assert header['projectScalarSchema']=='native-cjson-explicit-special-values-v6'
                source_records=[strict(line) for line in lines[1:]]
                with tp_raw.open() as handle:reference=list(csv.DictReader(handle,delimiter='\t'))
                columns=strict(mp.read_text());assert len(source_records)==len(reference)
                assert [r['iter'] for r in source_records]==list(range(len(source_records)))
                if native['exit_code']==0:assert len(source_records)==21
                finite=[]
                for iteration,(row,tsv) in enumerate(zip(source_records,reference)):
                    actual=next(emitted,None);assert actual is not None
                    assert int(tsv[columns['iter']])==iteration
                    value=row['parameters//']['S1/']['ASRV.Gamma:alpha']
                    token=tsv[columns['S1/ASRV.Gamma:alpha']]
                    native_value=float(token)
                    state='finite'
                    if type(value) is str:
                        assert value=='__project_scalar_v6__:positive_infinity'
                        assert native_value==math.inf
                        tags=row['numericParameterQuality//']['nonfinite']
                        assert tags==[dict(kind='positive_infinity',path=['S1/','ASRV.Gamma:alpha'])]
                        assert row['numericParameterQuality//']['literalNullPaths']==[]
                        value=None;state='positive_infinity';counts['positive_infinity_alpha_observations']+=1
                        summary['nonfinite_alpha_iterations'].append(iteration)
                    else:
                        assert type(value) in (int,float) and math.isfinite(value) and value>0
                        assert math.isclose(value,native_value,rel_tol=2e-13,abs_tol=0)
                        assert row['numericParameterQuality//']['nonfinite']==[]
                        finite.append(value)
                    expected=dict(stage=stage,chain_id=role['chain_id'],family=role['family'],
                        prior_label=role['prior_label'],iteration=iteration,alpha=value,alpha_state=state,
                        native_alpha_token=token,latent_log_alpha_observed=False)
                    for key in ['prior','likelihood','posterior']:
                        v=row['statistics//'][key]
                        assert type(v) in (int,float) and math.isfinite(v)
                        assert math.isclose(v,float(tsv[columns[key]]),rel_tol=2e-13,abs_tol=0)
                        expected[key]=v
                    actual_log=actual.pop('reconstructed_log_alpha')
                    assert actual==expected
                    if value is None:
                        assert actual_log is None;special.append({**expected,'reconstructed_log_alpha':None})
                    else:
                        independent_log=float(Decimal(str(value)).ln())
                        delta=abs(actual_log-independent_log);maximum_log_difference=max(maximum_log_difference,delta)
                        assert math.isclose(actual_log,independent_log,rel_tol=2e-13,abs_tol=2e-13)
                    counts['scalar_rows']+=1
                summary['rows']=len(source_records);summary['maximum_finite_alpha']=max(finite) if finite else None
                roles.append(summary)
        assert next(emitted,None) is None
    assert counts==producer['counts'] and counts['scalar_rows']==34020
    assert roles==producer['role_summaries'] and special==producer['special_observations']
    assert len(programs)==producer['continuous_generated_priors_checked']==411
    native=Path('results/ancestral/full-current-log-alpha-diagnostic-20261005-v2/deterministic-native-probe/attempt-0001')
    rows=[strict(line) for line in (native/'stdout.log').read_text().splitlines()]
    assert len(rows)==producer['native_deterministic_cases']==57
    grid=[-20,-10,-6,-3,0,6,20,50,100,500,700,709,709.7,709.78,709.79,710,750,1000,5000]
    overflows=0
    with localcontext() as context:
        context.prec=90
        maximum_double=Decimal(2)**1024-Decimal(2)**971
        for label,(mu,scale) in priors.items():
            cases=[row for row in rows if row['prior']==label]
            assert [row['log_alpha'] for row in cases]==grid
            for row in cases:
                x=Decimal(str(row['log_alpha']));alpha=x.exp()
                density=-(2*Decimal(scale)).ln()-abs(x-Decimal(mu))/Decimal(scale)
                assert row['mu']==mu and row['scale']==scale
                assert row['latent_log_density_tag']=='finite'
                assert math.isclose(row['latent_log_density'],float(density),rel_tol=2e-13,abs_tol=2e-13)
                if alpha>maximum_double:
                    assert row['alpha_tag']=='positive_infinity'
                    assert row['alpha']=='__log_alpha_diag_v2__:positive_infinity'
                    assert row['derived_log_density_tag']=='negative_infinity'
                    assert row['derived_log_density']=='__log_alpha_diag_v2__:negative_infinity'
                    overflows+=1
                else:
                    assert row['alpha_tag']=='finite' and type(row['alpha']) in (int,float)
                    assert abs(Decimal(str(row['alpha']))/alpha-1)<=Decimal('2e-13')
                    assert row['derived_log_density_tag']=='finite'
                    assert math.isclose(row['derived_log_density'],float(density-x),rel_tol=2e-13,abs_tol=2e-13)
                rates=row['gamma_rates']
                assert len(rates)==4 and all(type(rate) in (int,float) and math.isfinite(rate) and rate>=0 for rate in rates)
                mean=sum(Decimal(str(rate)) for rate in rates)/4
                assert abs(mean-1)<=Decimal('1e-10')
                if x>=100:assert rates==[1,1,1,1]
    assert overflows==producer['native_exp_overflow_cases']==15
    for path in [pp,tp,Path(__file__)]:pins[str(path)]=sha(path)
    for path,digest in pins.items():assert sha(path)==digest,path
    result=dict(status='passed_independent_full_current_alpha_trace_and_native_decimal_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(),scalar_rows=34020,role_slots=1644,
        retained_original_native_failures=24,positive_infinity_alpha_observations=26,
        review_roles=14,continuous_generated_priors_checked=411,native_deterministic_cases=57,
        native_overflow_cases=15,independent_decimal_precision_digits=90,
        maximum_reconstructed_log_alpha_difference=maximum_log_difference,
        producer_receipt_sha256=sha(pp),source_hashes=pins,scientific_eligibility=False,
        reviewed_arrays_admitted=False,original_chain_cause_proven=False,new_mcmc_runs=0,gpu=False,
        scope='Independent standard JSON/TSV parser reconstructs every available alpha and finite '
              'prior/likelihood/posterior field across all1644role slots, retains24original failures '
              'and checks full emitted rows/role summaries/review identities. Decimal90digit logarithms, '
              'exponentiation, finite-double boundary and analytic Laplace/pushforward densities verify '
              'all57native cases; rate normalization and large-shape/overflow unit behavior checked. '
              'General finite-shape Gamma category accuracy is not newly qualified. Original latent '
              'values unavailable: no historical causal proof, review admission or adequate posterior.')
    with args.receipt.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
