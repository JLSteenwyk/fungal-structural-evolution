"""Audit all qualification changes and coverage accounting without dropping draws."""
import copy
import json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha,write_json
from qualify_matched_kr_continuations import qualification
from calibrate_matched_kr import evaluate_refit,summarize
from matched_calibration_intervals import summarize_replicates
from audit_matched_kr_simulations import same
from plot_matched_kr_simulations import coverage_table


def main():
    threadpool_limits(1)
    root=Path('results/model_validation/matched-kr-continuation-overlay-20260928-v1')
    original=Path('results/model_validation/matched-kr-simulations-20260928-v1')
    oldaudit=Path('results/model_validation/matched-kr-simulation-audit-20260928-v1')
    r=json.loads((root/'receipt.json').read_text())
    assert r['source_receipt_sha256']==sha(original/'receipt.json')
    assert r['source_audit_receipt_sha256']==sha(oldaudit/'receipt.json')
    for name,digest in r['artifacts'].items():assert sha(root/name)==digest
    plan=json.loads((original/'run_plan.json').read_text())
    originals=json.loads((original/'manifest.json').read_text())
    overlays=json.loads((root/'overlay_manifest.json').read_text())
    mapped={(o['case'],o['replicate']):o for o in overlays}
    assert len(mapped)==len(overlays)==184
    assert set(mapped)=={(e['case'],e['replicate']) for e in originals if e['refit_status']=='refit_requires_review'}
    kr={};base={};strict_unmet=0
    for case in plan['cases']:
        with np.load(case['design'],allow_pickle=False) as h:a={k:h[k] for k in h.files}
        intervals=[];refits=[]
        for entry in [e for e in originals if e['case']==case['id']]:
            assert sha(entry['path'])==entry['sha256']
            saved=json.loads(Path(entry['path']).read_text());refit=saved['refit']
            key=(entry['case'],entry['replicate'])
            if key not in mapped:
                intervals.append(saved['interval']);refits.append(refit);continue
            e=mapped[key];assert sha(root/e['path'])==e['sha256']
            new=json.loads((root/e['path']).read_text())
            assert new['source_path']==entry['path'] and new['source_sha256']==entry['sha256']
            retries={}
            for path,digest in new['retry_bindings'].items():
                assert sha(path)==digest;retry=json.loads(Path(path).read_text());retries[retry['start']]=retry
            q=qualification(refit['fit'],retries);assert q==new['qualification'] and q['qualified']
            fixture_fit=refit['fit']
            strict_unmet+=not q['stricter_retry_threshold_all_met']
            adapter=dict(refit,status='refit_numerically_checked')
            recomputed=evaluate_refit(adapter,a['background'],a['family'],a['factor'],a['design'],a['beta'])
            same(recomputed,new['interval']);assert new['selected_parameters_unchanged']
            from scipy.stats import t
            f=refit['fit'];z=(np.array(f['beta'])-a['beta'])/np.sqrt(np.diag(f['conditional_beta_covariance']))
            covers=(abs(z)<=t.ppf(.975,f['residual_degrees_of_freedom'])).tolist()
            assert covers==new['conditional_t_covers_truth']
            adapter['nominal_095_t_interval_covers_truth']=covers
            intervals.append(recomputed);refits.append(adapter)
        assert len(intervals)==len(refits)==999
        kr[case['id']]=summarize(intervals);base[case['id']]=summarize_replicates(refits,case['coefficients'])
    same(kr,json.loads((root/'kr_coverage.json').read_text()))
    same(base,json.loads((root/'conditional_t_coverage.json').read_text()))
    table=coverage_table(plan,kr,base)
    old=coverage_table(plan,json.loads((oldaudit/'kr_coverage.json').read_text()),json.loads((oldaudit/'conditional_t_coverage.json').read_text()))
    keys=['case','method','coefficient'];paired=table.merge(old,on=keys,suffixes=('_new','_old'),validate='one_to_one')
    assert len(paired)==112 and (table.unresolved==0).all()
    assert (paired.observed_fraction_lower_new>=paired.observed_fraction_lower_old).all()
    assert (paired.observed_fraction_upper_new<=paired.observed_fraction_upper_old).all()
    assert (paired.lower_new>=paired.lower_old).all() and (paired.upper_new<=paired.upper_old).all()
    # Negative checks against a real qualification record.
    rejected=0
    for mutation in ['gradient','success','objective']:
        bad=copy.deepcopy(retries);endpoint=next(iter(bad.values()))['attempts'][-1]
        if mutation=='gradient':endpoint['projected_gradient_maximum']=.1
        elif mutation=='success':endpoint['optimizer_success']=False
        else:endpoint['objective']+=1
        assert not qualification(fixture_fit,bad)['qualified'];rejected+=1
    out=Path('results/model_validation/matched-kr-continuation-audit-20260928-v1');out.mkdir(parents=True,exist_ok=False)
    table.to_csv(out/'all_coefficient_coverage.tsv',sep='\t',index=False)
    result=dict(status='passed_all_184_qualification_overlays_and_15984_denominator_readback',
        attempted=15984,newly_qualified=184,unchanged_original_qualified=15800,
        coefficient_method_rows=112,unresolved_coefficient_outcomes=0,
        stricter_retry_diagnostic_unmet=strict_unmet,invalid_retry_evidence_rejected=rejected,
        original_coverage_bounds_contain_all_revised_results=True,
        source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='Numerical qualification update on the same frozen responses and parameters, '
              'using original tolerances. Original histories and stricter diagnostic preserved. '
              'No universal interval coverage, global optimum or biological claim.')
    write_json(out/'receipt.json',result);write_json(Path('metadata/matched_kr_continuation_audit_20260928.json'),result)
    print(result)


if __name__=='__main__':main()
