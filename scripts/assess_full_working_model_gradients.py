#!/usr/bin/env python3
"""Assess analytic stationarity across every completed production disposition."""
import argparse,hashlib,json,subprocess,time
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from cached_matched_likelihood import CachedMatchedLikelihood
from matched_reml_gradient import evaluate_gradient
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    config=json.loads(args.plan.read_text());ch=sha(args.plan)
    def verify():
        assert sha(args.plan)==ch
        for path,h in config['pins'].items():assert sha(path)==h,path
    verify();dependency=json.loads(Path(config['audit_launch']).read_text())
    while True:
        try:
            p=psutil.Process(dependency['pid'])
            if p.create_time()!=dependency['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dependency['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',dependency['unit'],'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines())
    assert state=={'ActiveState':'inactive','ExecMainStatus':'0'},state
    plan_path=Path(config['production_plan']);plan=json.loads(plan_path.read_text());ph=sha(plan_path)
    root=Path(plan['output']);audit=json.loads(Path(config['audit']).read_text());receipt=json.loads((root/'receipt.json').read_text())
    assert audit['status']=='passed_full_working_model_output_integrity_audit' and audit['source_receipt_sha256']==sha(root/'receipt.json')
    assert receipt['plan_sha256']==ph
    for name,h in receipt['artifacts'].items():assert sha(root/name)==h,name
    verify()
    entries={}
    for line in (root/'fit_manifest.jsonl').open():
        entry=json.loads(line);key=entry['fit_input_id'],entry['tree'];assert key not in entries;entries[key]=entry
    assert len(entries)==144040
    recipes=[json.loads(line) for line in (Path(plan['inventory'])/'unique_fit_recipes.jsonl').open()]
    nodes=pd.read_csv(plan['nodes'],sep='\t');targets=nodes[nodes.role.eq('target')][['node_id','guide','family_component']].rename(columns={'node_id':'target_id'})
    pairs=pd.read_csv(plan['pairs'],sep='\t',usecols=['target_id','background_id','species_pattern_id']).merge(targets,on='target_id',validate='many_to_one')
    pairs['row_identity']=[hashlib.sha256(json.dumps(list(r),separators=(',',':')).encode()).hexdigest() for r in pairs[['target_id','background_id','family_component','species_pattern_id']].itertuples(index=False,name=None)]
    selections=pd.read_csv(plan['selections'],sep='\t',usecols=['target_id','background_id','domain_config_id','policy','scenario_id']).merge(pairs,on=['target_id','background_id'],validate='many_to_one').sort_values('target_id',kind='stable')
    assert len(selections)==2786912
    pattern_index=pd.read_csv(Path(plan['factors'])/'patterns.tsv',sep='\t').set_index('species_pattern_id').row_index
    factors={p.stem:np.load(p)['factor'] for p in Path(plan['factors']).glob('*.npz')};assert len(factors)==5
    assert set(entries)=={(r['fit_input_id'],tree) for r in recipes for tree in factors}
    numeric=['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference']
    source=Path(plan['summaries']);out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    seen=set();counts=Counter();transitions=Counter();maximum=0.
    with (out/'analytic_stationarity.jsonl').open('w') as handle:
        for part in json.loads((source/'partition_manifest.json').read_text()):
            wanted={(r['guide'],r['policy'],r['scenario_id']):r for r in recipes if r['partition_index']==part['index']}
            if not wanted:continue
            path=source/part['path'];assert sha(path)==part['sha256']
            records=selections.merge(pd.read_parquet(path,columns=['domain_config_id']+numeric),on='domain_config_id',validate='many_to_one',sort=False)
            for group,frame in records.groupby(['guide','policy','scenario_id'],sort=True):
                if group not in wanted:continue
                recipe=wanted[group];identifier=recipe['fit_input_id']
                assert frame.target_id.is_monotonic_increasing and len(frame)==recipe['records']
                matrix=np.ascontiguousarray(frame[numeric].to_numpy(),dtype='<f8');matrix[matrix==0]=0.
                assert hashlib.sha256(matrix.tobytes()).hexdigest()==recipe['values_sha256']
                assert hashlib.sha256(np.asarray(frame.row_identity,dtype='S64').tobytes()).hexdigest()==recipe['ordered_identity_sha256']
                active=np.ptp(matrix[:,1:],axis=0)>1e-12;scales=np.std(matrix[:,1:][:,active],axis=0)
                x=np.column_stack([np.ones(len(matrix)),matrix[:,1:][:,active]/scales]);y=matrix[:,0]
                bg=pd.factorize(frame.background_id,sort=True)[0];fam=pd.factorize(frame.family_component,sort=True)[0]
                indices=frame.species_pattern_id.map(pattern_index);assert indices.notna().all();indices=indices.to_numpy(dtype=int)
                for tree,factor in factors.items():
                    key=identifier,tree;assert key not in seen;seen.add(key)
                    entry=entries[key];assert sha(entry['path'])==entry['sha256']
                    saved=json.loads(Path(entry['path']).read_text());assert saved['plan_sha256']==ph
                    payload=saved['payload'];assert payload['status']==entry['status']
                    result=dict(fit_input_id=identifier,tree=tree,original_status=payload['status'],source_fit_sha256=entry['sha256'])
                    if payload['status']=='fit_error_requires_review':
                        result.update(assessment='original_fit_error_unassessed',error_type=payload['error_type'])
                    else:
                        cache=CachedMatchedLikelihood(bg,fam,factor[indices],x,y)
                        theta=np.array(payload['log1p_ratios']);upper=np.log1p(payload['maximum_ratio'])
                        try:analytic=evaluate_gradient(cache,np.expm1(theta))
                        except (ArithmeticError,np.linalg.LinAlgError) as error:
                            result.update(assessment='analytic_diagnostic_unresolved',error_type=type(error).__name__,error=str(error))
                        else:
                            np.testing.assert_allclose(analytic['negative_profiled_reml'],payload['negative_profiled_reml'],rtol=1e-9,atol=1e-7)
                            maximum=max(maximum,abs(analytic['negative_profiled_reml']-payload['negative_profiled_reml']))
                            gradient=analytic['log1p_ratio_gradient']
                            projected=np.where(theta<=1e-7,np.minimum(gradient,0),np.where(theta>=upper-1e-7,np.maximum(gradient,0),gradient))
                            passed=bool(np.max(abs(projected))<=1e-3)
                            flags=payload['checks'];other=flags['best_optimizer_success'] and not flags['upper_bound_contact'] and flags['all_full_face_starts_agree']
                            result.update(assessment='analytic_score_evaluated',analytic_gradient=gradient.tolist(),analytic_projected_gradient=projected.tolist(),
                                          analytic_threshold_pass=passed,other_optimizer_checks_pass=bool(other),original_gradient_pass=flags['projected_gradient_pass'],
                                          all_numerical_checks_with_analytic_gradient=bool(passed and other))
                            transitions[str(flags['projected_gradient_pass'])+'->'+str(passed)]+=1
                    counts[result['assessment']]+=1;handle.write(json.dumps(result,allow_nan=False)+'\n')
                    if len(seen)%1000==0:handle.flush();print('Assessed full analytic stationarity',len(seen),'/ 144040',flush=True)
    assert seen==set(entries)
    verify()
    result=dict(status='complete_full_analytic_stationarity_assessment_pending_readback',plan_sha256=ch,
                dispositions=len(seen),assessment_counts=dict(counts),gradient_threshold_transitions=dict(transitions),maximum_objective_error=maximum,
                production_receipt_sha256=sha(root/'receipt.json'),production_audit_sha256=sha(config['audit']),
                artifacts={'analytic_stationarity.jsonl':sha(out/'analytic_stationarity.jsonl')},
                scope='All production dispositions retained and source inputs reverified. Original statuses unchanged. Same gradient threshold with independently fixture-validated analytic derivative; unresolved and original error cases explicit. Numerical diagnostics only, not global-optimum proof, calibrated inference or biological model adequacy.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
