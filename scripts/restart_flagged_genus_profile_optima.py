#!/usr/bin/env python3
"""Restart flagged unconstrained MG94 fits from all eight audited profile solutions."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time
from audit_genus_branch_parameter_profiles import declarations, unchanged_model_text, sha, table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['profiles','audit','plan','install','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    plan=json.loads(a.plan.read_text());audit=json.loads((a.audit/'receipt.json').read_text());config=json.loads((a.profiles/'config.json').read_text())
    if sha(a.profiles/'receipt.json')!=audit['source_receipt_sha256'] or sha(a.audit/'receipt.json')!=plan['audit_receipt_sha256']:raise ValueError('Audit provenance differs')
    for name,digest in audit['artifacts'].items():
        if sha(a.audit/name)!=digest:raise ValueError('Audited artifact changed')
    selected=[r for r in table(a.audit/'case_summary.tsv') if float(r['best_grid_exceeds_unconstrained_by'])>plan['discrepancy_threshold']]
    if sorted(r['case_id'] for r in selected)!=plan['case_ids']:raise ValueError('Flagged case universe differs')
    exe=a.install/'bin/hyphy'
    if sha(exe)!=config['executable_sha256']:raise ValueError('Executable differs')
    source_receipts={r['case_id']:r['receipt_sha256'] for r in json.loads((a.profiles/'receipt.json').read_text())['case_receipts']}
    a.output.mkdir(parents=True);rows=[];summaries=[];proofs=[]
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    for review in selected:
        case=review['case_id'];source=a.profiles/case
        if sha(source/'receipt.json')!=source_receipts[case]:raise ValueError('Case receipt changed')
        receipt=json.loads((source/'receipt.json').read_text());bylabel={r['fit_kind']:r for r in receipt['proofs']}
        original_unconstrained=next(r['log_likelihood'] for r in receipt['rows'] if r['fit_kind']=='unconstrained_reoptimization')
        previous_best=max(r['log_likelihood'] for r in receipt['rows']);case_rows=[]
        for start_row in receipt['rows']:
            factor=start_row['branch_parameter_multiplier'];label='unconstrained' if factor=='' else 'factor_'+str(factor).replace('.','p')
            src=source/label
            for name,digest in bylabel[label]['artifacts'].items():
                if sha(src/name)!=digest:raise ValueError('Saved starting solution changed')
            saved=(src/'profile_fit.bf').read_text();parameters,constrained=declarations(saved)
            names=[n for n in parameters if n.endswith('.tree_0.'+review['target_node']+'.t')]
            if len(names)!=1:raise ValueError('Ambiguous target')
            target=names[0]
            if constrained!=(set() if factor=='' else {target}):raise ValueError('Unexpected source constraints')
            start=saved.replace(target+':=',target+'=')
            values,fixed=declarations(start)
            if values!=parameters or fixed or unchanged_model_text(start,parameters)!=unchanged_model_text(saved,parameters):raise ValueError('Release changed source model')
            lfs=re.findall(r'^LikelihoodFunction\s+(\S+)\s*=',start,re.M)
            if len(lfs)!=1:raise ValueError('Ambiguous likelihood')
            lf=lfs[0];work=a.output/case/label;work.mkdir(parents=True)
            (work/'start.bf').write_text(start)
            code='ExecuteAFile('+json.dumps(str((work/'start.bf').resolve()))+');\nOPTIMIZATION_PRECISION=0.000001;\nUSE_LAST_RESULTS=1;\nOptimize(result,'+lf+');\nLFCompute('+lf+',LF_START_COMPUTE);\nLFCompute('+lf+',post_ll);\nLFCompute('+lf+',LF_DONE_COMPUTE);\n'
            code+='fprintf(stdout,"OPTIMUM\\t",Format(result[1][0],0,16),"\\t",Format(post_ll,0,16),"\\n");\n'
            for name in parameters:code+='fprintf(stdout,"PARAM\\t'+name+'\\t",Format('+name+',0,16),"\\n");\n'
            (work/'optimize.bf').write_text(code);began=time.monotonic()
            with (work/'optimize.log').open('w') as log:subprocess.run([str(exe.resolve()),'CPU=1','ENV=RANDOM_SEED=123;',str((work/'optimize.bf').resolve())],stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
            lines=(work/'optimize.log').read_text().splitlines();opt=[x.split('\t')[1:] for x in lines if x.startswith('OPTIMUM\t')];pairs=[x.split('\t')[1:] for x in lines if x.startswith('PARAM\t')];fitted={n:float(v) for n,v in pairs}
            if len(opt)!=1 or len(fitted)!=len(pairs) or set(fitted)!=set(parameters) or any(not math.isfinite(v) or v<0 for v in fitted.values()):raise ValueError('Invalid optimized parameter grid')
            ll,post=map(float,opt[0])
            if not math.isfinite(ll) or abs(ll-post)>1e-6 or ll<start_row['log_likelihood']-1e-5:raise ValueError('Invalid or worsened optimum')
            exported=start
            for name,value in fitted.items():
                pattern=r'^(global )?'+re.escape(name)+r'=[-+0-9.eE]+;'
                replacement=('global ' if '.model_MGREV.' in name else '')+name+'='+repr(value)+';'
                exported,n=re.subn(pattern,lambda m:replacement,exported,flags=re.M)
                if n!=1:raise ValueError('Ambiguous export')
            replay_values,replay_fixed=declarations(exported)
            if replay_values!=fitted or replay_fixed or unchanged_model_text(exported,parameters)!=unchanged_model_text(saved,parameters):raise ValueError('Export changed model content')
            (work/'fit.bf').write_text(exported)
            code='ExecuteAFile('+json.dumps(str((work/'fit.bf').resolve()))+');\nLFCompute('+lf+',LF_START_COMPUTE);\nLFCompute('+lf+',read_ll);\nLFCompute('+lf+',LF_DONE_COMPUTE);\nfprintf(stdout,"READBACK=",Format(read_ll,0,16),"\\n");\n'
            (work/'readback.bf').write_text(code)
            with (work/'readback.log').open('w') as log:subprocess.run([str(exe.resolve()),'CPU=1',str((work/'readback.bf').resolve())],stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
            replay=[float(x.partition('=')[2]) for x in (work/'readback.log').read_text().splitlines() if x.startswith('READBACK=')]
            if len(replay)!=1 or abs(replay[0]-ll)>1e-6:raise ValueError('Fresh model likelihood differs')
            row={'case_id':case,'start_solution':label,'starting_log_likelihood':start_row['log_likelihood'],'log_likelihood':ll,'gain_over_start':ll-start_row['log_likelihood'],'gain_over_previous_unconstrained':ll-original_unconstrained,'gain_over_previous_best_grid_or_unconstrained':ll-previous_best,'target_t':fitted[target],'omega':next(v for n,v in fitted.items() if n.endswith('.omega')),'fresh_readback_error':replay[0]-ll,'elapsed_seconds':time.monotonic()-began}
            rows.append(row);case_rows.append(row);proofs.append({'case_id':case,'start':label,'artifacts':{str(x.relative_to(a.output)):sha(x) for x in work.iterdir()}})
        best=max(case_rows,key=lambda r:r['log_likelihood'])
        summaries.append({'case_id':case,'starts':len(case_rows),'best_start':best['start_solution'],'best_log_likelihood':best['log_likelihood'],'best_target_t':best['target_t'],'best_omega':best['omega'],'gain_over_previous_unconstrained':best['gain_over_previous_unconstrained'],'gain_over_previous_best_evaluated':best['gain_over_previous_best_grid_or_unconstrained'],'range_of_restarted_log_likelihoods':max(r['log_likelihood'] for r in case_rows)-min(r['log_likelihood'] for r in case_rows),'fcs_exposed_case':review['fcs_exposed_case']})
        print('Restarted',case,flush=True)
    if len(rows)!=plan['fits']:raise ValueError('Incomplete restart grid')
    for name,data in [('restart_points.tsv',rows),('case_summary.tsv',summaries)]:
        with (a.output/name).open('w',newline='') as handle:
            writer=csv.DictWriter(handle,list(data[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(data)
    result={'status':'complete_flagged_unconstrained_multistart_diagnostics','cases':len(summaries),'fits':len(rows),'fresh_model_readbacks':len(rows),'maximum_readback_error':max(abs(r['fresh_readback_error']) for r in rows),'source_profile_receipt_sha256':sha(a.profiles/'receipt.json'),'source_audit_receipt_sha256':sha(a.audit/'receipt.json'),'plan_sha256':sha(a.plan),'script_sha256':sha(Path(__file__)),'audit_helper_sha256':sha(Path(__file__).with_name('audit_genus_branch_parameter_profiles.py')),'proofs':proofs,'artifacts':{name:sha(a.output/name) for name in ['restart_points.tsv','case_summary.tsv']},'interpretation':'Unconstrained restarts from all eight previously evaluated solutions for four cases flagged by grid/unconstrained discrepancy. Original outputs retained. Better evaluated solutions do not prove global optimality. Fixed-grid profiles have not been regenerated from these new optima; no dS intervals, saturation clearance or selection inference.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()
