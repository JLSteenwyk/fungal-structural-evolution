#!/usr/bin/env python3
"""Verify complete native-rate replay provenance and compare numerical errors."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from ancestral_chain_attempt import sha,write_json


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show','fungal-fastml-native-rate-replay-20260928.service','-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),state
    producer=json.loads(Path(plan['producer_plan']).read_text());checked=set()
    def check(path,h):
        assert sha(path)==h,path;checked.add(str(path))
    for config in [plan,producer]+[j['config'] for j in producer['jobs'] if j['config'] is not None]:
        for path,h in config['pins'].items():check(path,h)
    root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['plan_sha256']==ph and receipt['fits']==306
    prior_root=Path(plan['readback']);prior=json.loads((prior_root/'receipt.json').read_text())
    rate_root=Path(plan['rates']);rate_receipt=json.loads((rate_root/'receipt.json').read_text())
    for name,h in rate_receipt['artifacts'].items():check(rate_root/name,h)
    rates={r['id']:r for r in map(json.loads,(rate_root/'rate_comparison.jsonl').read_text().splitlines())}
    jobs={j['id']:j for j in producer['jobs']};expected={k for k,j in jobs.items() if j['job']['character_count']}
    assert set(receipt['artifacts'])=={key+'.json' for key in expected}
    assert set(receipt['empty_inputs'])==set(jobs)-expected and len(receipt['empty_inputs'])==6
    rows=[]
    for name,h in receipt['artifacts'].items():
        check(root/name,h);r=json.loads((root/name).read_text());key=r['id'];assert name==key+'.json'
        entry=prior['readbacks'][key];check(entry['path'],entry['sha256']);old=json.loads(Path(entry['path']).read_text())
        assert r['source_readback_sha256']==entry['sha256'] and r['job_id']==old['job_id'] and r['variant']==old['variant']
        assert r['probability_rows']==old['probability_rows'] and not r['raw_values_clipped']
        assert r['reported_log_likelihood']==old['reported_log_likelihood']
        assert r['scipy_rate_likelihood_difference']==old['likelihood_difference']
        assert r['scipy_rate_maximum_probability_difference']==old['maximum_probability_difference']
        assert r['native_rates']==rates[key]['native_rates']
        assert r['native_rate_likelihood_difference']==r['native_rate_replay_log_likelihood']-r['reported_log_likelihood']
        native=json.loads((Path(producer['output'])/key/'readback.json').read_text());rp=Path(native['attempt_receipt'])
        check(rp,r['source_attempt_sha256']);assert r['source_attempt_sha256']==native['attempt_receipt_sha256']
        for artifact,digest in json.loads(rp.read_text())['artifacts'].items():check(rp.parent/artifact,digest)
        for k,v in r.items():
            if isinstance(v,float):assert np.isfinite(v),(key,k)
        rows.append(r)
    table=pd.DataFrame(rows);assert len(table)==len(set(table.id))==306
    summaries={}
    for variant,frame in table.groupby('variant'):
        assert len(frame)==153 and int(frame.probability_rows.sum())==8058340
        summaries[variant]=dict(fits=len(frame),probability_rows=int(frame.probability_rows.sum()),
            scipy_rate_maximum_absolute_likelihood_difference=float(frame.scipy_rate_likelihood_difference.abs().max()),
            native_rate_maximum_absolute_likelihood_difference=float(frame.native_rate_likelihood_difference.abs().max()),
            scipy_rate_maximum_probability_difference=float(frame.scipy_rate_maximum_probability_difference.max()),
            native_rate_maximum_probability_difference=float(frame.native_rate_maximum_probability_difference.max()),
            native_likelihood_difference_over_1e_minus7=int((frame.native_rate_likelihood_difference.abs()>1e-7).sum()))
    paired=[]
    for job_id,frame in table.groupby('job_id'):
        assert len(frame)==2
        a=frame.set_index('variant').loc['precision-only'];b=frame.set_index('variant').loc['precision-cache-refresh']
        paired.append(dict(job_id=job_id,native_replay_likelihood_change=float(b.native_rate_replay_log_likelihood-a.native_rate_replay_log_likelihood),
                           precision_only_report_error=float(a.native_rate_likelihood_difference),cache_refresh_report_error=float(b.native_rate_likelihood_difference)))
    paired=pd.DataFrame(paired);worse=paired[paired.native_replay_likelihood_change < -1e-7]
    args.output.mkdir(parents=True,exist_ok=False)
    table.drop(columns='native_rates').to_csv(args.output/'replay_comparison.tsv',sep='\t',index=False)
    paired.to_csv(args.output/'paired_native_likelihoods.tsv',sep='\t',index=False)
    result=dict(status='complete_native_rate_replay_provenance_and_scalar_readback',source_receipt_sha256=sha(root/'receipt.json'),
        variants=summaries,checked_files=len(checked),cache_refresh_worse_than_precision_only=worse.to_dict('records'),
        remaining_precision_only_report_discrepancies=table.loc[table.variant.eq('precision-only') & table.native_rate_likelihood_difference.abs().gt(1e-7),['job_id','native_rate_likelihood_difference']].to_dict('records'),
        pins={str(p):sha(p) for p in [args.plan,Path(__file__)]},artifacts={p.name:sha(p) for p in args.output.iterdir()},
        scope='Verified full expected fit set, original input/configuration pins, native attempt artifacts, rate identity and scalar arithmetic. Uses completed independent-pruning replay rather than another inference implementation. 1e-7 is descriptive discrepancy reporting, not proof of global optimization or biological adequacy.')
    write_json(args.output/'receipt.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
