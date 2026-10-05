#!/usr/bin/env python3
"""Close full alpha diagnostic/reader/figure transports and replay figure bins."""
import argparse
from collections import defaultdict
import csv
from datetime import datetime,timezone
import json
import math
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();assert not args.receipt.exists()
    pins={};closed=[]
    stages=[('ancestral_log_alpha_probe_debugger','v1',79425),
            ('ancestral_log_alpha_diagnostic','v2',38416),
            ('ancestral_log_alpha_readback','v1',92677),
            ('ancestral_log_alpha_figure','v1',82726)]
    for prefix,version,session in stages:
        vp=Path(f'metadata/{prefix}_20261005_{version}.json')
        ep=Path(f'metadata/{prefix}_execution_20261005_{version}.json')
        tp=Path(f'metadata/{prefix}_transport_20261005_{version}.json')
        pp=Path(f'metadata/{prefix}_original_tool_payloads_20261005_{version}.json')
        v,e,t,tool=[json.loads(path.read_text()) for path in [vp,ep,tp,pp]]
        assert e['status']=='exited_zero_with_receipt' and e['exit_code']==0 and not e['timed_out']
        assert e['receipt_sha256']==t['validation_sha256']==sha(vp)
        assert t['original_tool_session_id']==tool['original_tool_session_id']==session
        assert tool['terminal']['exit_code']==0
        assert t['whole_wrapper_initial_and_terminal_payloads_matched']
        assert t['manager_start_records']==t['manager_completion_records']==1
        assert v['scientific_eligibility'] is False
        for mapping in [e['source_hashes'],e['artifacts'],t['source_hashes'],v['source_hashes']]:
            for path,digest in mapping.items():bind(pins,path,digest)
        for path in [vp,ep,tp,pp]:bind(pins,path)
        closed.append(dict(stage=prefix,version=version,original_tool_session_id=session,
                           wrapper_exit_code=0,whole_original_execution_verified=True))
    diagnostic=json.loads(Path('metadata/ancestral_log_alpha_diagnostic_20261005_v2.json').read_text())
    reader=json.loads(Path('metadata/ancestral_log_alpha_readback_20261005_v1.json').read_text())
    assert reader['scalar_rows']==diagnostic['counts']['scalar_rows']==34020
    assert reader['positive_infinity_alpha_observations']==diagnostic['counts']['positive_infinity_alpha_observations']==26
    assert reader['role_slots']==1644 and reader['review_roles']==14
    assert reader['native_deterministic_cases']==57 and reader['native_overflow_cases']==15
    failure_path=Path('metadata/ancestral_log_alpha_diagnostic_failure_audit_20261005_v1.json')
    failure=json.loads(failure_path.read_text())
    assert failure['original_tool_terminal_exit_code']==1 and failure['deterministic_native_exit_code']==-11
    assert failure['native_output_rows']==0 and failure['whole_original_wrapper_payloads_matched']
    for path,digest in failure['source_hashes'].items():bind(pins,path,digest)
    bind(pins,failure_path)
    visual_path=Path('metadata/ancestral_log_alpha_figure_visual_review_20261005_v1.json')
    visual=json.loads(visual_path.read_text())
    for path,digest in visual['source_hashes'].items():bind(pins,path,digest)
    bind(pins,visual_path)
    backtrace=Path('metadata/ancestral_log_alpha_probe_original_backtrace_20261005_v1.log')
    assert sha(backtrace)==sha('results/ancestral/log-alpha-probe-debugger-20261005-v1/attempt-0001/stdout.log')
    assert 'builtin_function_ejson_null' in backtrace.read_text()
    bind(pins,backtrace)
    figure=json.loads(Path('metadata/ancestral_log_alpha_figure_20261005_v1.json').read_text())
    bins=defaultdict(list)
    with Path('results/ancestral/full-current-log-alpha-diagnostic-20261005-v2/all_current_alpha_traces.jsonl').open() as handle:
        for line in handle:
            row=json.loads(line);bins[(row['stage'],row['prior_label'],row['iteration'])].append(row)
    with Path(figure['table']).open() as handle:table=list(csv.DictReader(handle,delimiter='\t'))
    assert len(table)==len(bins)==126
    maximum=0.
    seen=set()
    for row in table:
        key=(row['stage'],row['prior'],int(row['iteration']))
        assert key not in seen;seen.add(key);items=bins[key]
        finite=sorted(math.log10(item['alpha']) for item in items if item['alpha_state']=='finite')
        assert int(row['role_rows'])==len(items)==(532 if key[0]=='v6' else 8)
        assert int(row['finite_rows'])==len(finite)
        assert int(row['positive_infinity_rows'])==len(items)-len(finite)
        for field,p in [('finite_log10_alpha_q10',.1),('finite_log10_alpha_median',.5),('finite_log10_alpha_q90',.9)]:
            position=p*(len(finite)-1);low=math.floor(position);high=math.ceil(position)
            value=finite[low]+(position-low)*(finite[high]-finite[low])
            actual=float(row[field]);maximum=max(maximum,abs(actual-value))
            assert math.isclose(actual,value,rel_tol=2e-13,abs_tol=2e-13)
    assert seen==set(bins)
    bind(pins,Path(__file__));verify(pins)
    result=dict(status='complete_verified_full_current_alpha_diagnostic_and_figure',
        checked_utc=datetime.now(timezone.utc).isoformat(),role_slots=1644,scalar_rows=34020,
        retained_original_native_failures=24,positive_infinity_alpha_observations=26,review_roles=14,
        generated_prior_programs=411,deterministic_native_cases=57,native_exp_overflow_cases=15,
        figure_bins_independently_replayed=126,maximum_figure_quantile_difference=maximum,
        original_successful_stage_transports=closed,first_probe_failure_retained=True,
        complete_bound_files=len(pins),source_hashes=pins,scientific_eligibility=False,
        original_latent_log_alpha_observed=False,original_chain_cause_proven=False,
        scalar_review_arrays_admitted=False,adequate_ancestral_posterior=False,
        all_eight_aims_incomplete=True,new_mcmc_runs=0,gpu=False,
        scope='Full current V6/V7comparison alpha trace/native prior readback and90digitDecimal '
              'native overflow/density controls close with actual original producer/reader/debugger/figure '
              'tool terminals and whole wrapper/manager journals.126figure bins manually sorted/replayed; '
              'public image copies match visually reviewed output. First native JSONnull conversion fault '
              'and failed wrapper remain preserved. Demonstrated compatible overflow mechanism does not '
              'recover historical latent values, prove chain cause, repair installed libraries, qualify '
              'general Gamma category rules, admit reviews or establish an adequate posterior/biological aim.')
    with args.receipt.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
