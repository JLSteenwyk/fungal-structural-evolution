#!/usr/bin/env python3
"""Reoptimize MG94 nuisance parameters across fixed longest-branch parameter grids."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def table(p):
    with p.open() as f: return list(csv.DictReader(f, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['fits','slices','plan','install','output']: parser.add_argument('--'+key,type=Path,required=True)
    a=parser.parse_args(); plan=json.loads(a.plan.read_text()); script_sha=sha(Path(__file__))
    for folder,key in [(a.fits,'fit_receipt_sha256'),(a.slices,'slice_receipt_sha256')]:
        if sha(folder/'receipt.json')!=plan[key]: raise ValueError('Changed source receipt')
    fits=json.loads((a.fits/'receipt.json').read_text()); slices=json.loads((a.slices/'receipt.json').read_text())
    for name,digest in slices['artifacts'].items():
        if sha(a.slices/name)!=digest: raise ValueError('Changed slice artifact')
    reviews=table(a.slices/'case_review.tsv'); source_cases={r['case_id']:r['receipt_sha256'] for r in fits['case_receipts']}
    fixed={(r['case_id'],float(r['branch_multiplier'])):float(r['log_likelihood']) for r in table(a.slices/'likelihood_slices.tsv')}
    if len(reviews)!=plan['cases'] or len(source_cases)!=len(reviews) or {r['case_id'] for r in reviews}!=set(source_cases): raise ValueError('Wrong case universe')
    exe=a.install/'bin/hyphy'
    if sha(exe)!=plan['executable_sha256']: raise ValueError('Changed executable')
    if a.output.exists(): raise FileExistsError('Use a new immutable output')
    a.output.mkdir(parents=True)
    cfg={'plan_sha256':sha(a.plan),'script_sha256':script_sha,'source_fit_receipt_sha256':sha(a.fits/'receipt.json'),'source_slice_receipt_sha256':sha(a.slices/'receipt.json'),'executable_sha256':sha(exe)}
    (a.output/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    def run(review):
        case=review['case_id'];src=a.fits/case;rc=json.loads((src/'receipt.json').read_text())
        if sha(src/'receipt.json')!=source_cases[case]: raise ValueError('Changed fit receipt')
        for name,digest in rc['artifacts'].items():
            if sha(src/name)!=digest: raise ValueError('Changed fit artifact')
        saved=(src/'fit.bf').read_text();lfs=re.findall(r'^LikelihoodFunction\s+(\S+)\s*=',saved,re.M)
        declarations=re.findall(r'^(global )?([A-Za-z0-9_.]+)=([-+0-9.eE]+);',saved,re.M)
        parameters={name:float(value) for prefix,name,value in declarations if (prefix and '.model_MGREV.' in name) or '.tree_0.' in name and name.endswith('.t')}
        targets=[name for name in parameters if name.endswith('.tree_0.'+review['target_node']+'.t')]
        if len(lfs)!=1 or len(targets)!=1 or len([n for n in parameters if '.model_MGREV.' in n])!=6: raise ValueError('Unexpected saved parameterization')
        expected_nodes=set(json.loads((src/'fit.json').read_text())['branch attributes']['0'])
        observed_nodes={n.split('.tree_0.')[1][:-2] for n in parameters if '.tree_0.' in n}
        if observed_nodes!=expected_nodes: raise ValueError('Saved parameter grid omits or adds branches')
        lf=lfs[0];target=targets[0];original=parameters[target];out=a.output/case;out.mkdir();records=[];proofs=[]
        for factor in [None]+plan['branch_parameter_multipliers']:
            label='unconstrained' if factor is None else 'factor_'+str(factor).replace('.','p');work=out/label;work.mkdir()
            bf='ExecuteAFile('+json.dumps(str((src/'fit.bf').resolve()))+');\n'
            if factor is not None: bf+=target+':='+repr(original*factor)+';\n'
            bf+='OPTIMIZATION_PRECISION='+str(plan['optimization_precision'])+';\nUSE_LAST_RESULTS=1;\nOptimize(profile_result,'+lf+');\nLFCompute('+lf+',LF_START_COMPUTE);\nLFCompute('+lf+',post_ll);\nLFCompute('+lf+',LF_DONE_COMPUTE);\n'
            bf+='fprintf(stdout,"OPTIMUM\\t",Format(profile_result[1][0],0,16),"\\t",Format(post_ll,0,16),"\\n");\n'
            for name in parameters: bf+='fprintf(stdout,"PARAM\\t'+name+'\\t",Format('+name+',0,16),"\\n");\n'
            (work/'optimize.bf').write_text(bf);start=time.monotonic()
            with (work/'optimize.log').open('w') as f: subprocess.run([str(exe.resolve()),'CPU=1','ENV=RANDOM_SEED=123;',str((work/'optimize.bf').resolve())],stdout=f,stderr=subprocess.STDOUT,env=env,check=True)
            log=(work/'optimize.log').read_text();lls=re.findall(r'^OPTIMUM\t(\S+)\t(\S+)',log,re.M);vals=re.findall(r'^PARAM\t(\S+)\t(\S+)',log,re.M)
            if len(lls)!=1 or len(vals)!=len(parameters): raise ValueError('Incomplete optimization output')
            values={name:float(value) for name,value in vals};ll,post=map(float,lls[0]);baseline=rc['log_likelihood'] if factor is None else fixed[case,float(factor)]
            if set(values)!=set(parameters) or any(not math.isfinite(v) or v<0 for v in values.values()) or not math.isfinite(ll) or abs(ll-post)>1e-6: raise ValueError('Invalid optimized parameters/likelihood')
            if factor is not None and not math.isclose(values[target],original*factor,rel_tol=1e-10,abs_tol=1e-10): raise ValueError('Fixed target moved')
            if ll<baseline-1e-5: raise ValueError('Optimizer reduced likelihood below its starting slice')
            # Export the entire model with all fitted free values, then reload it in a new process.
            exported=saved
            for name,value in values.items():
                pattern=r'^(global )?'+re.escape(name)+r'=[-+0-9.eE]+;'
                replacement=(('global ' if '.model_MGREV.' in name else '')+name+(':=' if name==target and factor is not None else '=')+repr(value)+';')
                exported,n=re.subn(pattern,lambda m:replacement,exported,flags=re.M)
                if n!=1: raise ValueError('Ambiguous fitted parameter export')
            (work/'profile_fit.bf').write_text(exported)
            replay='ExecuteAFile('+json.dumps(str((work/'profile_fit.bf').resolve()))+');\nLFCompute('+lf+',LF_START_COMPUTE);\nLFCompute('+lf+',readback_ll);\nLFCompute('+lf+',LF_DONE_COMPUTE);\nfprintf(stdout,"READBACK=",Format(readback_ll,0,16),"\\n");\n'
            (work/'readback.bf').write_text(replay)
            with (work/'readback.log').open('w') as f: subprocess.run([str(exe.resolve()),'CPU=1',str((work/'readback.bf').resolve())],stdout=f,stderr=subprocess.STDOUT,env=env,check=True)
            matches=re.findall(r'^READBACK=(\S+)',(work/'readback.log').read_text(),re.M)
            if len(matches)!=1 or abs(float(matches[0])-ll)>1e-6: raise ValueError('Fresh saved-fit likelihood differs')
            row={'case_id':case,'target_node':review['target_node'],'fit_kind':'unconstrained_reoptimization' if factor is None else 'fixed_branch_parameter_profile','branch_parameter_multiplier':'' if factor is None else factor,'target_parameter':values[target],'log_likelihood':ll,'delta_ll_from_original_fit':ll-rc['log_likelihood'],'improvement_over_fixed_nuisance_slice':ll-baseline,'omega':next(v for n,v in values.items() if n.endswith('.omega')),'changed_nuisance_parameters':sum(not math.isclose(values[n],v,rel_tol=1e-6,abs_tol=1e-9) for n,v in parameters.items() if n!=target),'fresh_readback_ll_error':float(matches[0])-ll,'elapsed_seconds':time.monotonic()-start,'original_target_dS':review['target_dS'],'marker_copy_caveat':review['marker_copy_caveat']}
            records.append(row);proofs.append({'fit_kind':label,'artifacts':{p.name:sha(p) for p in work.iterdir() if p.is_file()}})
        best=max(r['log_likelihood'] for r in records)
        for row in records: row['delta_ll_from_best_evaluated_fit']=row['log_likelihood']-best
        receipt={'status':'complete_branch_parameter_profile_with_fresh_fit_readbacks','case_id':case,'source_fit_receipt_sha256':sha(src/'receipt.json'),'rows':records,'proofs':proofs}
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        return records,{'case_id':case,'receipt_sha256':sha(out/'receipt.json')}
    rows=[];proofs=[]
    pool=ThreadPoolExecutor(max_workers=plan['workers'])
    futures=[pool.submit(run,r) for r in reviews]
    try:
        for i,future in enumerate(as_completed(futures),1):
            data,proof=future.result();rows.extend(data);proofs.append(proof)
            if i%25==0: print('Profiled',i,'of',len(reviews),flush=True)
    except BaseException:
        for future in futures: future.cancel()
        raise
    finally:
        pool.shutdown(wait=True,cancel_futures=True)
    rows.sort(key=lambda r:(r['case_id'],r['fit_kind'],str(r['branch_parameter_multiplier'])))
    with (a.output/'profile_points.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    if sha(Path(__file__))!=script_sha or sha(a.plan)!=cfg['plan_sha256']: raise ValueError('Producer or plan changed during execution')
    r={'status':'complete_full_case_branch_parameter_profiles_pending_summary_audit','cases':len(proofs),'optimized_fits':len(rows),'fresh_saved_fit_readbacks':len(rows),'config_sha256':sha(a.output/'config.json'),'case_receipts':sorted(proofs,key=lambda r:r['case_id']),'artifacts':{'profile_points.tsv':sha(a.output/'profile_points.tsv')},'interpretation':plan['interpretation']}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print('Completed',len(proofs),'cases',len(rows),'fits',flush=True)


if __name__=='__main__':main()
