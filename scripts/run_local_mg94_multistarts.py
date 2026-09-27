#!/usr/bin/env python3
"""Unconstrained MG94 refits from every audited profile solution, all cases."""
import argparse,csv,json,math,os,re,subprocess,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from audit_local_branch_parameter_profiles import declarations,unchanged_model_text
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify();source=Path(plan['starts']);sr=json.loads((source/'receipt.json').read_text())
    assert sr['cases']==1632 and sr['starts']==13056
    for name,h in sr['artifacts'].items():assert sha(source/name)==h
    rows=list(csv.DictReader((source/'start_manifest.tsv').open(),delimiter='\t'))
    seeds={(s['case_id'],s['start_label']):s for s in map(json.loads,(source/'parameter_seeds.jsonl').open())}
    groups={}
    for row in rows:groups.setdefault(row['case_id'],[]).append(row)
    assert len(groups)==1632 and all(len(g)==8 for g in groups.values()) and len(seeds)==13056
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    (out/'run_plan.json').write_bytes(args.plan.read_bytes())
    executable=Path(plan['executable']).resolve();env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    def invoke(script,log):
        with log.open('w') as handle:
            subprocess.run([str(executable),'CPU=1','ENV=RANDOM_SEED=123;',str(script.resolve())],stdout=handle,stderr=subprocess.STDOUT,env=env,check=True)
    def run(case,starts):
        case_out=out/case;case_out.mkdir();results=[];proofs=[]
        for row in starts:
            label=row['start_label'];work=case_out/label;work.mkdir()
            original_path=Path(row['original_free_model']);seed_path=Path(row['profile_seed_model'])
            assert sha(original_path)==row['original_free_model_sha256'] and sha(seed_path)==row['profile_seed_model_sha256']
            original=original_path.read_text();parameters,fixed=declarations(original);assert not fixed
            seed=seeds[case,label];values=seed['parameters'];saved_seed,saved_fixed=declarations(seed_path.read_text())
            assert values==saved_seed and set(values)==set(parameters) and len(values)==int(row['parameters'])
            assert unchanged_model_text(original,parameters)==unchanged_model_text(seed_path.read_text(),parameters)
            target=row['target_parameter_name'];assert target in values
            assert saved_fixed==(set() if label=='unconstrained' else {target})
            functions=re.findall(r'^LikelihoodFunction\s+(\S+)\s*=',original,re.M);assert len(functions)==1;lf=functions[0]
            # Start from the free model; assigning values never imports profile := constraints.
            code='ExecuteAFile('+json.dumps(str(original_path.resolve()))+');\n'
            code+=''.join(name+'='+repr(value)+';\n' for name,value in sorted(values.items()))
            code+='LFCompute('+lf+',LF_START_COMPUTE);\nLFCompute('+lf+',initial_ll);\nLFCompute('+lf+',LF_DONE_COMPUTE);\nfprintf(stdout,"INITIAL\\t",Format(initial_ll,0,16),"\\n");\n'
            code+='OPTIMIZATION_PRECISION='+str(plan['optimization_precision'])+';\nUSE_LAST_RESULTS=1;\nOptimize(refit_result,'+lf+');\nLFCompute('+lf+',LF_START_COMPUTE);\nLFCompute('+lf+',post_ll);\nLFCompute('+lf+',LF_DONE_COMPUTE);\n'
            code+='fprintf(stdout,"OPTIMUM\\t",Format(refit_result[1][0],0,16),"\\t",Format(post_ll,0,16),"\\n");\n'
            code+=''.join('fprintf(stdout,"PARAM\\t'+name+'\\t",Format('+name+',0,16),"\\n");\n' for name in sorted(values))
            (work/'optimize.bf').write_text(code);started=time.monotonic();invoke(work/'optimize.bf',work/'optimize.log')
            log=(work/'optimize.log').read_text();initial=re.findall(r'^INITIAL\t(\S+)',log,re.M);opt=re.findall(r'^OPTIMUM\t(\S+)\t(\S+)',log,re.M);pars=re.findall(r'^PARAM\t(\S+)\t(\S+)',log,re.M)
            assert len(initial)==len(opt)==1 and len(pars)==len(values)
            fitted={name:float(value) for name,value in pars};ll,post=map(float,opt[0]);start_ll=float(initial[0])
            assert set(fitted)==set(values) and all(math.isfinite(v) and v>=0 for v in fitted.values())
            assert math.isfinite(ll) and abs(ll-post)<=1e-6 and abs(start_ll-float(row['starting_log_likelihood']))<=1e-6
            assert ll>=start_ll-1e-5
            exported=original
            for name,value in fitted.items():
                pattern=r'^(global )?'+re.escape(name)+r'=[-+0-9.eE]+;'
                replacement=('global ' if '.model_MGREV.' in name else '')+name+'='+repr(value)+';'
                exported,count=re.subn(pattern,lambda match:replacement,exported,flags=re.M);assert count==1
            exported_values,exported_fixed=declarations(exported)
            assert exported_values==fitted and not exported_fixed
            assert unchanged_model_text(original,parameters)==unchanged_model_text(exported,parameters)
            (work/'fit.bf').write_text(exported)
            replay='ExecuteAFile('+json.dumps(str((work/'fit.bf').resolve()))+');\nLFCompute('+lf+',LF_START_COMPUTE);\nLFCompute('+lf+',check_ll);\nLFCompute('+lf+',LF_DONE_COMPUTE);\nfprintf(stdout,"READBACK=",Format(check_ll,0,16),"\\n");\n'
            (work/'readback.bf').write_text(replay);invoke(work/'readback.bf',work/'readback.log')
            matches=re.findall(r'^READBACK=(\S+)',(work/'readback.log').read_text(),re.M);assert len(matches)==1
            readback=float(matches[0]);assert math.isfinite(readback) and abs(readback-ll)<=1e-6
            results.append(dict(case_id=case,start_label=label,starting_log_likelihood=start_ll,log_likelihood=ll,gain_over_start=ll-start_ll,
                                target_parameter=fitted[target],omega=next(v for n,v in fitted.items() if n.endswith('.omega')),
                                changed_parameters=sum(not math.isclose(fitted[n],v,rel_tol=1e-6,abs_tol=1e-9) for n,v in values.items()),
                                readback_error=readback-ll,elapsed_seconds=time.monotonic()-started))
            proofs.append(dict(start_label=label,source_original_sha256=row['original_free_model_sha256'],source_seed_sha256=row['profile_seed_model_sha256'],
                               parameter_seed_sha256=row['parameter_sha256'],artifacts={p.name:sha(p) for p in work.iterdir()}))
        best=max(r['log_likelihood'] for r in results)
        for row in results:row['deficit_from_best_start']=best-row['log_likelihood']
        result=dict(status='complete_all_eight_unconstrained_starts_with_fresh_readbacks',case_id=case,rows=results,proofs=proofs)
        (case_out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
        return results,dict(case_id=case,receipt_sha256=sha(case_out/'receipt.json'))
    all_rows=[];case_receipts=[]
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
        futures=[pool.submit(run,case,starts) for case,starts in sorted(groups.items())]
        try:
            for ix,future in enumerate(as_completed(futures),1):
                data,proof=future.result();all_rows.extend(data);case_receipts.append(proof)
                if ix%10==0:print('Completed all unconstrained starts',ix,'/ 1632',flush=True)
        except BaseException:
            for future in futures:future.cancel()
            raise
    assert len(all_rows)==13056 and len(case_receipts)==1632
    all_rows.sort(key=lambda r:(r['case_id'],r['start_label']))
    with (out/'fit_summary.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,list(all_rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(all_rows)
    verify()
    result=dict(status='complete_full_unconstrained_multistarts_pending_independent_audit',cases=1632,fits=13056,plan_sha256=ph,
                case_receipts=sorted(case_receipts,key=lambda r:r['case_id']),artifacts={name:sha(out/name) for name in ['run_plan.json','fit_summary.tsv']},
                scope='Every case and all eight audited profile seeds, all model parameters free. Initial, optimized and freshly reloaded likelihoods checked. Not global-optimum proof, branch dS confidence intervals, selection eligibility or biological conclusion.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
