#!/usr/bin/env python3
"""Restartable complete five-tree working-model grid, gated by full source audit."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import fcntl
import hashlib
import json
import multiprocessing
from pathlib import Path
import subprocess
import time
import traceback
import numpy as np
import pandas as pd
import psutil
from threadpoolctl import threadpool_limits
from screen_duplication_domain_alignment_coverage import sha
from fit_matched_mixed_model_cached import fit_variance_ratios
from matched_mixed_covariance import MatchedCovariance, profiled_reml

FACTORS={}


def initialize(factor_root):
    global FACTORS
    threadpool_limits(limits=1)
    for path in sorted(Path(factor_root).glob('*.npz')):
        with np.load(path) as arrays:FACTORS[path.stem]=arrays['factor']
    assert len(FACTORS)==5


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def fit_case(task):
    identifier, bg, family, pattern_rows, matrix, output, plan_hash=task
    directory=Path(output)/'fits'/identifier[:2]/identifier
    directory.mkdir(parents=True,exist_ok=True)
    outcomes=[]
    for tree,full_factor in FACTORS.items():
        path=directory/(tree+'.json')
        if path.exists():
            saved=json.loads(path.read_text())
            assert saved['plan_sha256']==plan_hash and saved['fit_input_id']==identifier and saved['tree']==tree
            assert saved['payload_sha256']==digest(saved['payload'])
            outcomes.append(dict(fit_input_id=identifier,tree=tree,status=saved['payload']['status'],path=str(path),sha256=sha(path)))
            continue
        try:
            covariates=matrix[:,1:]
            active=np.ptp(covariates,axis=0)>1e-12
            assert np.all(abs(covariates[:,~active])<=1e-12), 'Nonzero constant covariate changes zero-contrast interpretation'
            scales=np.std(covariates[:,active],axis=0)
            x=np.column_stack([np.ones(len(matrix)),covariates[:,active]/scales])
            y=matrix[:,0]
            factor=full_factor[pattern_rows]
            fitted=fit_variance_ratios(bg,family,factor,x,y)
            # Re-evaluate all candidate parameters using the separately checked direct evaluator.
            maximum_error=0.; direct_cache={}
            for candidate in fitted['candidates']:
                key=tuple(candidate['theta'])
                if key not in direct_cache:
                    ratios=np.expm1(candidate['theta'])
                    direct_cache[key]=profiled_reml(MatchedCovariance(bg,family,factor,1.,*ratios),x,y)
                reference=direct_cache[key]['negative_profiled_reml']
                np.testing.assert_allclose(reference,candidate['objective'],rtol=1e-9,atol=1e-7)
                maximum_error=max(maximum_error,abs(reference-candidate['objective']))
            best=profiled_reml(MatchedCovariance(bg,family,factor,1.,*fitted['ratios']),x,y)
            for name in ['beta','profiled_scale','residual_quadratic','conditional_beta_covariance']:
                np.testing.assert_allclose(fitted[name],best[name],rtol=1e-7,atol=1e-8)
            conversion=np.concatenate([[1.],1./scales])
            fitted.update(records=len(matrix),active_covariates=active.tolist(),covariate_scales=scales.tolist(),
                          raw_unit_beta=(np.asarray(fitted['beta'])*conversion).tolist(),
                          raw_unit_conditional_beta_covariance=(np.asarray(fitted['conditional_beta_covariance'])*conversion[:,None]*conversion[None,:]).tolist(),
                          direct_candidate_readback=dict(candidates_checked=len(fitted['candidates']),distinct_parameter_vectors=len(direct_cache),maximum_objective_error=float(maximum_error)),
                          scope='Conditional working-model candidate; no calibrated confidence interval, significance test or causal duplication claim. All original settings remain mapped externally.')
            payload=fitted
        except Exception as error:
            payload=dict(status='fit_error_requires_review',error_type=type(error).__name__,error=str(error),traceback=traceback.format_exc())
        saved=dict(plan_sha256=plan_hash,fit_input_id=identifier,tree=tree,payload=payload,payload_sha256=digest(payload))
        temporary=path.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(saved,indent=2,allow_nan=False)+'\n')
        temporary.replace(path)
        outcomes.append(dict(fit_input_id=identifier,tree=tree,status=payload['status'],path=str(path),sha256=sha(path)))
    return outcomes


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify()
    dependency=json.loads(Path(plan['inventory_audit_launch']).read_text())
    while True:
        try:
            process=psutil.Process(dependency['pid'])
            if process.create_time()!=dependency['created'] or process.status()==psutil.STATUS_ZOMBIE:break
            assert process.cmdline()==dependency['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',dependency['unit'],'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines())
    assert state=={'ActiveState':'inactive','ExecMainStatus':'0'},state
    audit=json.loads(Path(plan['inventory_audit']).read_text())
    inventory=Path(plan['inventory']);ir=json.loads((inventory/'receipt.json').read_text())
    assert audit['status']=='passed_full_matched_fit_inventory_readback' and audit['source_receipt_sha256']==sha(inventory/'receipt.json')
    assert audit['unique_record_inputs']==ir['unique_record_inputs']==28808
    verify()
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True)
    lock=(out/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    saved_plan=out/'run_plan.json'
    if saved_plan.exists():assert saved_plan.read_bytes()==args.plan.read_bytes()
    else:saved_plan.write_bytes(args.plan.read_bytes())
    binding=out/'source_audit_binding.json'
    binding_text=json.dumps(dict(inventory_audit_sha256=sha(plan['inventory_audit']),inventory_receipt_sha256=sha(inventory/'receipt.json')),indent=2)+'\n'
    if binding.exists():assert binding.read_text()==binding_text
    else:binding.write_text(binding_text)
    recipes=[json.loads(line) for line in (inventory/'unique_fit_recipes.jsonl').open()]
    by_setting={(r['partition_index'],r['guide'],r['policy'],r['scenario_id']):r for r in recipes}
    assert len(by_setting)==len(recipes)
    nodes=pd.read_csv(plan['nodes'],sep='\t')
    targets=nodes[nodes.role.eq('target')][['node_id','guide','family_component']].rename(columns={'node_id':'target_id'})
    pairs=pd.read_csv(plan['pairs'],sep='\t',usecols=['target_id','background_id','species_pattern_id']).merge(targets,on='target_id',validate='many_to_one')
    pairs['row_identity']=[hashlib.sha256(json.dumps(list(row),separators=(',',':')).encode()).hexdigest() for row in pairs[['target_id','background_id','family_component','species_pattern_id']].itertuples(index=False,name=None)]
    selections=pd.read_csv(plan['selections'],sep='\t',usecols=['target_id','background_id','policy','scenario_id','domain_config_id'])
    selections=selections.merge(pairs,on=['target_id','background_id'],validate='many_to_one').sort_values('target_id',kind='stable')
    assert len(selections)==2786912
    pattern_table=pd.read_csv(Path(plan['factors'])/'patterns.tsv',sep='\t').set_index('species_pattern_id')
    source=Path(plan['summaries']);partitions=json.loads((source/'partition_manifest.json').read_text())
    numeric=['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference']
    pending=set();submitted=set();counts=Counter();completed=0;error_count=0
    with (out/'fit_manifest.jsonl').open('w') as manifest, ProcessPoolExecutor(max_workers=plan['workers'],mp_context=multiprocessing.get_context('spawn'),initializer=initialize,initargs=(plan['factors'],)) as pool:
        def collect(block=False):
            nonlocal pending,completed,error_count
            done,pending=wait(pending,return_when=FIRST_COMPLETED,timeout=None if block else 0)
            for future in done:
                for result in future.result():
                    counts[result['status']]+=1;completed+=1
                    error_count+=result['status']=='fit_error_requires_review'
                    manifest.write(json.dumps(result,sort_keys=True)+'\n')
                manifest.flush()
                print('Completed working-model dispositions',completed,'/ 144040',dict(counts),flush=True)
            if error_count>=10 and error_count/max(1,completed)>.1:
                raise RuntimeError('More than ten percent numerical errors after at least ten errors; outputs preserved for review')
        for part in partitions:
            needed={key[1:]:recipe for key,recipe in by_setting.items() if key[0]==part['index']}
            if not needed:continue
            path=source/part['path'];assert sha(path)==part['sha256']
            values=pd.read_parquet(path,columns=['domain_config_id']+numeric)
            records=selections.merge(values,on='domain_config_id',validate='many_to_one',sort=False)
            for key,frame in records.groupby(['guide','policy','scenario_id'],sort=True):
                if key not in needed:continue
                recipe=needed[key];identifier=recipe['fit_input_id']
                assert identifier not in submitted and frame.target_id.is_monotonic_increasing
                matrix=np.ascontiguousarray(frame[numeric].to_numpy(),dtype='<f8');matrix[matrix==0]=0.
                assert len(frame)==recipe['records'] and np.isfinite(matrix).all()
                assert hashlib.sha256(matrix.tobytes()).hexdigest()==recipe['values_sha256']
                assert hashlib.sha256(np.asarray(frame.row_identity,dtype='S64').tobytes()).hexdigest()==recipe['ordered_identity_sha256']
                bg=pd.factorize(frame.background_id,sort=True)[0];family=pd.factorize(frame.family_component,sort=True)[0]
                pattern_rows=frame.species_pattern_id.map(pattern_table.row_index)
                assert pattern_rows.notna().all()
                pending.add(pool.submit(fit_case,(identifier,bg,family,pattern_rows.to_numpy(dtype=int),matrix,str(out),ph)))
                submitted.add(identifier)
                if len(pending)>=2*plan['workers']:collect(block=True)
            collect()
        while pending:collect(block=True)
    assert submitted=={r['fit_input_id'] for r in recipes} and completed==144040
    verify()
    receipt=dict(status='complete_full_working_model_dispositions_pending_full_audit',plan_sha256=ph,
                 unique_inputs=len(submitted),tree_fit_dispositions=completed,status_counts=dict(counts),
                 artifacts={p.name:sha(p) for p in [out/'run_plan.json',out/'source_audit_binding.json',out/'fit_manifest.jsonl']},
                 scope='All exact unique inputs and five trees attempted with all 22 starts/faces and direct candidate likelihood readback. Review/error cases retained. Original 82944 settings map externally. Not inferential calibration, model adequacy, causal evidence or project completion.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
