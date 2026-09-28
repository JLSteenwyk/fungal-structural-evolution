#!/usr/bin/env python3
"""Refine every flagged working-model fit without altering production outputs."""
import argparse, fcntl, hashlib, json, time, traceback
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from refine_matched_reml_analytic import refine
from screen_duplication_domain_alignment_coverage import sha


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    config=json.loads(args.plan.read_text());ch=sha(args.plan)
    threadpool_limits(1)
    plan_path=Path(config['production_plan']);plan=json.loads(plan_path.read_text());ph=sha(plan_path)
    def verify():
        assert sha(args.plan)==ch and sha(plan_path)==ph
        for path,h in {**plan['pins'],**config['pins']}.items():assert sha(path)==h,path
    verify()
    review_root=Path(config['review_root'])
    receipt=json.loads((review_root/'receipt.json').read_text())
    assert receipt['status']=='passed_full_analytic_stationarity_output_readback'
    cases=[json.loads(line) for line in (review_root/'remaining_review_cases.jsonl').read_text().splitlines()]
    entries={(r['fit_input_id'],r['tree']):r for r in cases}
    assert len(entries)==len(cases)==config['expected_cases']
    wanted_ids={r['fit_input_id'] for r in cases}
    recipes=[json.loads(line) for line in (Path(plan['inventory'])/'unique_fit_recipes.jsonl').open()]
    recipes=[r for r in recipes if r['fit_input_id'] in wanted_ids]
    assert {r['fit_input_id'] for r in recipes}==wanted_ids
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=True)
    lock=(out/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    snapshot=out/'run_plan.json'
    if snapshot.exists():assert snapshot.read_bytes()==args.plan.read_bytes()
    else:snapshot.write_bytes(args.plan.read_bytes())
    nodes=pd.read_csv(plan['nodes'],sep='\t');targets=nodes[nodes.role.eq('target')][['node_id','guide','family_component']].rename(columns={'node_id':'target_id'})
    pairs=pd.read_csv(plan['pairs'],sep='\t',usecols=['target_id','background_id','species_pattern_id']).merge(targets,on='target_id',validate='many_to_one')
    pairs['row_identity']=[hashlib.sha256(json.dumps(list(r),separators=(',',':')).encode()).hexdigest() for r in pairs[['target_id','background_id','family_component','species_pattern_id']].itertuples(index=False,name=None)]
    selections=pd.read_csv(plan['selections'],sep='\t',usecols=['target_id','background_id','domain_config_id','policy','scenario_id']).merge(pairs,on=['target_id','background_id'],validate='many_to_one').sort_values('target_id',kind='stable')
    assert len(selections)==2786912
    pattern_index=pd.read_csv(Path(plan['factors'])/'patterns.tsv',sep='\t').set_index('species_pattern_id').row_index
    factors={p.stem:np.load(p)['factor'] for p in Path(plan['factors']).glob('*.npz')};assert len(factors)==5
    numeric=['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference']
    source=Path(plan['summaries']);seen=set();counts=Counter();manifest=[];started=time.perf_counter()
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
                key=identifier,tree
                if key not in entries:continue
                assert key not in seen;seen.add(key)
                entry=entries[key];assert sha(entry['source_fit'])==entry['source_fit_sha256']
                saved=json.loads(Path(entry['source_fit']).read_text());assert saved['plan_sha256']==ph
                old=saved['payload']
                path=out/(identifier+'-'+tree+'.json')
                if path.exists():
                    result=json.loads(path.read_text())
                    assert result['plan_sha256']==ch and result['source']==entry
                    assert result['payload_sha256']==digest(result['payload'])
                else:
                    try:
                        assert np.all(abs(matrix[:,1:][:,~active])<=1e-12)
                        payload=refine(bg,fam,factor[indices],x,y,old['log1p_ratios'],old['maximum_ratio'],config['maxiter'])
                        np.testing.assert_allclose(payload['original_objective'],old['negative_profiled_reml'],rtol=1e-9,atol=1e-7)
                        conversion=np.r_[1.,1./scales]
                        payload.update(records=len(matrix),active_covariates=active.tolist(),covariate_scales=scales.tolist(),
                            raw_unit_beta=(np.asarray(payload['beta'])*conversion).tolist(),
                            raw_unit_conditional_beta_covariance=(np.asarray(payload['conditional_beta_covariance'])*conversion[:,None]*conversion[None,:]).tolist())
                    except Exception as error:
                        payload=dict(status='refinement_error_requires_review',error=str(error),error_type=type(error).__name__,traceback=traceback.format_exc())
                    result=dict(plan_sha256=ch,source=entry,payload=payload,payload_sha256=digest(payload))
                    temporary=path.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');temporary.replace(path)
                counts[result['payload']['status']]+=1
                manifest.append(dict(fit_input_id=identifier,tree=tree,path=str(path),sha256=sha(path),status=result['payload']['status']))
                print('Refined',len(seen),'/',len(entries),dict(counts),flush=True)
    assert seen==set(entries)
    verify()
    (out/'manifest.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in manifest))
    result=dict(status='complete_review_fit_refinement_pending_independent_readback',plan_sha256=ch,
                dispositions=len(seen),status_counts=dict(counts),wall_seconds=time.perf_counter()-started,
                artifacts={p.name:sha(p) for p in [out/'manifest.jsonl',out/'run_plan.json']},
                scope='All remaining flagged unique tree fits retained. Original grid unchanged. Numerical refinement requires independent output readback; biological adequacy and inferential calibration unresolved.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
