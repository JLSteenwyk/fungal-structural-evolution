"""Describe residual summaries in all existing 15984 Gaussian simulations."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha,write_json
from simulate_matched_working_model import simulate
from refit_matched_simulation import response_seed
from matched_marginal_residual_diagnostics import diagnostics


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    config=json.loads(args.plan.read_text());ch=sha(args.plan);threadpool_limits(1)
    def verify():
        assert sha(args.plan)==ch
        for name,h in config['pins'].items():assert sha(name)==h,name
    verify();plan=json.loads(Path(config['simulation_plan']).read_text())
    for name,h in plan['pins'].items():assert sha(name)==h,name
    root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    for name,h in receipt['artifacts'].items():assert sha(root/name)==h
    overlay=Path(config['overlay']);orr=json.loads((overlay/'receipt.json').read_text())
    audit=json.loads(Path(config['qualification_audit']).read_text())
    assert audit['status']=='passed_all_184_qualification_overlays_and_15984_denominator_readback'
    assert audit['source_receipt_sha256']==sha(overlay/'receipt.json')
    for name,h in orr['artifacts'].items():assert sha(overlay/name)==h
    qualified={}
    for entry in json.loads((overlay/'overlay_manifest.json').read_text()):
        path=overlay/entry['path'];assert sha(path)==entry['sha256']
        r=json.loads(path.read_text());assert r['qualification']['qualified'] and r['selected_parameters_unchanged']
        key=(entry['case'],entry['replicate']);assert key not in qualified;qualified[key]=r
    manifest=json.loads((root/'manifest.json').read_text())
    expected={(c['id'],rep) for c in plan['cases'] for rep in range(plan['replicates'])}
    assert len(manifest)==len(expected)==15984 and {(r['case'],r['replicate']) for r in manifest}==expected
    records=[];seen=set();overlay_used=set()
    for case in plan['cases']:
        with np.load(case['design'],allow_pickle=False) as h:a={k:h[k] for k in h.files}
        for entry in [r for r in manifest if r['case']==case['id']]:
            path=Path(entry['path']);assert sha(path)==entry['sha256']
            saved=json.loads(path.read_text());refit=saved['refit'];rep=entry['replicate'];key=(case['id'],rep)
            assert saved['case']==case and saved['replicate']==rep and refit['fit_id']==case['id']
            assert saved['plan_sha256']==sha(config['simulation_plan'])
            if refit['status']=='refit_requires_review':
                assert qualified[key]['source_sha256']==entry['sha256'];overlay_used.add(key)
            else:assert refit['status']=='refit_numerically_checked'
            entropy=response_seed(plan['master_seed'],case['id'],rep);assert refit['seed_entropy']==entropy
            y=simulate(a['background'],a['family'],a['factor'],a['design'],a['beta'],case['scale'],case['ratios'],1,np.random.default_rng(np.random.SeedSequence(entropy)))[:,0]
            assert hashlib.sha256(np.asarray(y,dtype='<f8').tobytes()).hexdigest()==refit['response_sha256']
            f=refit['fit']
            for mode,v in [('generating_covariance',case['scale']*np.r_[1.,case['ratios']]),('fitted_covariance',f['profiled_scale']*np.r_[1.,f['ratios']])]:
                result=diagnostics(a['background'],a['family'],a['factor'],a['design'],y,v)
                row=dict(case=case['id'],replicate=rep,mode=mode,status=result['status'],source_sha256=entry['sha256'],original_refit_status=refit['status'],qualified_by_continuation=key in qualified)
                if result['status']=='descriptive_marginal_residual_diagnostics':
                    if mode=='fitted_covariance':np.testing.assert_allclose(result['beta'],f['beta'],rtol=1e-7,atol=1e-8)
                    s=result['summaries'];row.update({k:v for k,v in s.items() if not isinstance(v,list)})
                    row['fraction_absolute_above_2']=s['absolute_above_2']/case['records']
                    row['fraction_absolute_above_3']=s['absolute_above_3']/case['records']
                    row['maximum_absolute_bin_mean']=max(abs(b['mean_standardized_residual']) for b in result['covariate_bins'])
                    row['maximum_bin_second_moment']=max(b['mean_squared_standardized_residual'] for b in result['covariate_bins'])
                else:row['reason']=result['reason']
                records.append(row)
            assert key not in seen;seen.add(key)
        print('Residual simulation case completed',case['id'],len(seen),'/15984',flush=True)
    assert seen==expected and overlay_used==set(qualified) and len(records)==31968
    frame=pd.DataFrame(records);summaries=[]
    metrics=['mean','second_raw_moment','third_raw_moment','fourth_raw_moment','fraction_absolute_above_2','fraction_absolute_above_3','maximum_absolute_bin_mean','maximum_bin_second_moment']
    for (case,mode),group in frame.groupby(['case','mode']):
        assert len(group)==999
        good=group[group.status=='descriptive_marginal_residual_diagnostics']
        for metric in metrics:
            values=good[metric]
            q=values.quantile([.025,.5,.975]) if len(values) else pd.Series([None]*3,index=[.025,.5,.975])
            summaries.append(dict(case=case,mode=mode,metric=metric,attempted=len(group),available=len(good),unresolved=len(group)-len(good),mean=values.mean(),q025=q.loc[.025],median=q.loc[.5],q975=q.loc[.975]))
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    frame.to_parquet(out/'replicate_diagnostics.parquet',index=False)
    pd.DataFrame(summaries).to_csv(out/'diagnostic_reference_summary.tsv',sep='\t',index=False)
    verify();write_json(out/'receipt.json',dict(status='complete_existing_simulation_residual_summaries_pending_readback',plan_sha256=ch,responses=len(seen),diagnostic_dispositions=len(frame),status_counts=frame.status.value_counts().to_dict(),
        artifacts={p.name:sha(p) for p in out.iterdir()},scope='Same15984 seeded responses and selected fits; no new simulations or refits. Generating and fitted covariance comparisons in two synthetic designs only. Empirical quantiles are descriptive, not transferable critical values or evidence of fungal model adequacy.'))

if __name__=='__main__':main()
