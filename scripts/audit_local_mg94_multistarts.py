#!/usr/bin/env python3
"""Audit all unconstrained codon starts and retain between-start discrepancies."""
import argparse,csv,json,math,re,subprocess,time
from pathlib import Path
import psutil
from audit_local_branch_parameter_profiles import declarations,unchanged_model_text
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser()
    for name in ['plan','launch','output']:ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    launch=json.loads(args.launch.read_text());lh=sha(args.launch);assert launch['plan_sha256']==ph
    while True:
        try:
            p=psutil.Process(launch['pid'])
            if p.create_time()!=launch['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines())
    assert state=={'ActiveState':'inactive','ExecMainStatus':'0'},state
    assert sha(args.plan)==ph and sha(args.launch)==lh
    for path,h in plan['pins'].items():assert sha(path)==h,path
    root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text());assert receipt['plan_sha256']==ph
    assert receipt['status']=='complete_full_unconstrained_multistarts_pending_independent_audit'
    for name,h in receipt['artifacts'].items():assert sha(root/name)==h,name
    source=Path(plan['starts'])
    starts={(r['case_id'],r['start_label']):r for r in csv.DictReader((source/'start_manifest.tsv').open(),delimiter='\t')}
    seeds={(r['case_id'],r['start_label']):r for r in map(json.loads,(source/'parameter_seeds.jsonl').open())}
    aggregate={(r['case_id'],r['start_label']):r for r in csv.DictReader((root/'fit_summary.tsv').open(),delimiter='\t')}
    assert set(starts)==set(seeds)==set(aggregate) and len(starts)==13056
    profiles=Path('results/cds/local-mg94-branch-parameter-profiles-20260927-v1')
    source_profile=json.loads((profiles/'receipt.json').read_text());source_cases={r['case_id']:r['receipt_sha256'] for r in source_profile['case_receipts']}
    summary=[];seen=set();hashes=parameters=0;maximum_error=0.
    for case_proof in receipt['case_receipts']:
        case=case_proof['case_id'];folder=root/case;assert sha(folder/'receipt.json')==case_proof['receipt_sha256']
        cr=json.loads((folder/'receipt.json').read_text());assert cr['case_id']==case
        assert cr['status']=='complete_all_eight_unconstrained_starts_with_fresh_readbacks' and len(cr['rows'])==len(cr['proofs'])==8
        proof_by_label={p['start_label']:p for p in cr['proofs']};assert len(proof_by_label)==8
        assert sha(profiles/case/'receipt.json')==source_cases[case]
        prior=json.loads((profiles/case/'receipt.json').read_text())['rows']
        best=max(r['log_likelihood'] for r in cr['rows']);worst=min(r['log_likelihood'] for r in cr['rows'])
        for row in cr['rows']:
            key=case,row['start_label'];assert key in starts and key not in seen;seen.add(key)
            start=starts[key];proof=proof_by_label[key[1]];work=folder/key[1]
            for name,h in proof['artifacts'].items():assert sha(work/name)==h;hashes+=1
            assert sha(start['original_free_model'])==start['original_free_model_sha256']==proof['source_original_sha256']
            assert sha(start['profile_seed_model'])==start['profile_seed_model_sha256']==proof['source_seed_sha256']
            original=Path(start['original_free_model']).read_text();original_parameters,original_fixed=declarations(original);assert not original_fixed
            seed_values,_=declarations(Path(start['profile_seed_model']).read_text());assert seed_values==seeds[key]['parameters']
            assert start['parameter_sha256']==proof['parameter_seed_sha256']==seeds[key]['parameter_sha256']
            code=(work/'optimize.bf').read_text()
            assert code.startswith('ExecuteAFile('+json.dumps(str(Path(start['original_free_model']).resolve()))+');\n')
            initial_values={n:float(v) for n,v in re.findall(r'^([A-Za-z0-9_.]+)=([-+0-9.eE]+);',code,re.M) if n in original_parameters}
            assert initial_values==seed_values and ':=' not in code
            exported=(work/'fit.bf').read_text();values,fixed=declarations(exported)
            assert not fixed and set(values)==set(original_parameters) and all(math.isfinite(v) and v>=0 for v in values.values())
            assert unchanged_model_text(exported,values)==unchanged_model_text(original,original_parameters)
            log=(work/'optimize.log').read_text();initial=re.findall(r'^INITIAL\t(\S+)',log,re.M);opt=re.findall(r'^OPTIMUM\t(\S+)\t(\S+)',log,re.M);pars=re.findall(r'^PARAM\t(\S+)\t(\S+)',log,re.M)
            assert len(initial)==len(opt)==1 and len(pars)==len(values)
            assert {n:float(v) for n,v in pars}==values;parameters+=len(values)
            replay=re.findall(r'^READBACK=(\S+)',(work/'readback.log').read_text(),re.M);assert len(replay)==1
            ll=row['log_likelihood'];start_ll=float(initial[0]);fresh=float(replay[0])
            assert row['starting_log_likelihood']==start_ll
            assert all(math.isfinite(v) for v in [ll,start_ll,fresh])
            assert abs(start_ll-float(start['starting_log_likelihood']))<=1e-6 and ll>=start_ll-1e-5
            assert max(abs(float(v)-ll) for v in opt[0])<=1e-6 and abs(fresh-ll)<=1e-6
            assert abs(row['readback_error']-(fresh-ll))<=1e-10 and abs(row['gain_over_start']-(ll-start_ll))<=1e-9
            assert abs(row['deficit_from_best_start']-(best-ll))<=1e-9
            assert row['target_parameter']==values[start['target_parameter_name']]
            assert row['omega']==next(v for n,v in values.items() if n.endswith('.omega'))
            assert row['changed_parameters']==sum(not math.isclose(values[n],v,rel_tol=1e-6,abs_tol=1e-9) for n,v in seed_values.items())
            assert all(str(v)==aggregate[key][k] for k,v in row.items())
            maximum_error=max(maximum_error,abs(fresh-ll))
        prior_unconstrained=next(r['log_likelihood'] for r in prior if r['fit_kind']=='unconstrained_reoptimization')
        prior_best=max(r['log_likelihood'] for r in prior)
        assert best>=prior_best-1e-5
        summary.append(dict(case_id=case,best_log_likelihood=best,worst_log_likelihood=worst,between_start_spread=best-worst,
                            starts_within_1e_5_of_best=sum(best-r['log_likelihood']<=1e-5 for r in cr['rows']),
                            best_gain_over_prior_unconstrained=best-prior_unconstrained,best_gain_over_prior_evaluated_fit=best-prior_best))
    assert seen==set(starts) and len(summary)==1632
    args.output.mkdir(parents=True,exist_ok=False)
    with (args.output/'case_summary.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,list(summary[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(sorted(summary,key=lambda r:r['case_id']))
    for path,h in plan['pins'].items():assert sha(path)==h,path
    result=dict(status='passed_full_local_mg94_multistart_artifact_and_numeric_audit',cases=1632,fits=13056,
                artifact_hashes_checked=hashes,parameters_checked=parameters,maximum_fresh_likelihood_error=maximum_error,
                cases_start_spread_gt_1e_5=sum(r['between_start_spread']>1e-5 for r in summary),
                source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),artifacts={'case_summary.tsv':sha(args.output/'case_summary.tsv')},
                scope='Full start identities, free constraints, parameters, initial/optimized/fresh saved likelihoods and case summaries checked. No fresh optimization rerun, proof of global optima, dS intervals or selection conclusion.')
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
