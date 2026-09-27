#!/usr/bin/env python3
"""Compare analytic and finite-difference scores for every fit in a frozen prefix."""
import hashlib,json
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from cached_matched_likelihood import CachedMatchedLikelihood
from matched_reml_gradient import evaluate_gradient
from screen_duplication_domain_alignment_coverage import sha


def main():
    snapshot_path=Path('metadata/full_working_model_early_diagnostics_20260927.json')
    snapshot=json.loads(snapshot_path.read_text())
    plan_path=Path('metadata/full_matched_working_model_plan_20260927.json');plan=json.loads(plan_path.read_text());ph=sha(plan_path)
    for path,h in plan['pins'].items():assert sha(path)==h,path
    proof_path=Path('metadata/matched_reml_gradient_checks_20260927.json');proof=json.loads(proof_path.read_text())
    assert proof['status']=='passed_analytic_matched_reml_gradient_dense_score_checks'
    for path,h in proof['pins'].items():assert sha(path)==h,path
    with (Path(plan['output'])/'fit_manifest.jsonl').open('rb') as handle:prefix=handle.read(snapshot['manifest_prefix_bytes'])
    assert hashlib.sha256(prefix).hexdigest()==snapshot['manifest_prefix_sha256']
    entries=[json.loads(line) for line in prefix.splitlines()];assert len(entries)==snapshot['completed']==320
    recipes={r['fit_input_id']:r for r in map(json.loads,(Path(plan['inventory'])/'unique_fit_recipes.jsonl').open())}
    selected_ids={e['fit_input_id'] for e in entries}
    nodes=pd.read_csv(plan['nodes'],sep='\t');targets=nodes[nodes.role.eq('target')][['node_id','guide','family_component']].rename(columns={'node_id':'target_id'})
    pairs=pd.read_csv(plan['pairs'],sep='\t',usecols=['target_id','background_id','species_pattern_id']).merge(targets,on='target_id',validate='many_to_one')
    pairs['row_identity']=[hashlib.sha256(json.dumps(list(r),separators=(',',':')).encode()).hexdigest() for r in pairs[['target_id','background_id','family_component','species_pattern_id']].itertuples(index=False,name=None)]
    selections=pd.read_csv(plan['selections'],sep='\t',usecols=['target_id','background_id','domain_config_id','policy','scenario_id']).merge(pairs,on=['target_id','background_id'],validate='many_to_one').sort_values('target_id',kind='stable')
    pattern_index=pd.read_csv(Path(plan['factors'])/'patterns.tsv',sep='\t').set_index('species_pattern_id').row_index
    factors={p.stem:np.load(p)['factor'] for p in Path(plan['factors']).glob('*.npz')}
    numeric=['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference']
    summaries=Path(plan['summaries']);prepared={}
    for part in json.loads((summaries/'partition_manifest.json').read_text()):
        wanted={(r['guide'],r['policy'],r['scenario_id']):r for k,r in recipes.items() if k in selected_ids and r['partition_index']==part['index']}
        if not wanted:continue
        path=summaries/part['path'];assert sha(path)==part['sha256']
        records=selections.merge(pd.read_parquet(path,columns=['domain_config_id']+numeric),on='domain_config_id',validate='many_to_one',sort=False)
        for key,frame in records.groupby(['guide','policy','scenario_id']):
            if key not in wanted:continue
            recipe=wanted[key];matrix=np.ascontiguousarray(frame[numeric].to_numpy(),dtype='<f8');matrix[matrix==0]=0.
            assert hashlib.sha256(matrix.tobytes()).hexdigest()==recipe['values_sha256']
            assert hashlib.sha256(np.asarray(frame.row_identity,dtype='S64').tobytes()).hexdigest()==recipe['ordered_identity_sha256']
            active=np.ptp(matrix[:,1:],axis=0)>1e-12;scales=np.std(matrix[:,1:][:,active],axis=0)
            x=np.column_stack([np.ones(len(matrix)),matrix[:,1:][:,active]/scales])
            prepared[recipe['fit_input_id']]=(pd.factorize(frame.background_id,sort=True)[0],pd.factorize(frame.family_component,sort=True)[0],frame.species_pattern_id.map(pattern_index).to_numpy(dtype=int),x,matrix[:,0])
    assert set(prepared)==selected_ids
    out=Path('results/model_validation/working-model-gradient-diagnostic-20260927-v1');out.mkdir(parents=True,exist_ok=False)
    counts=Counter();max_objective_error=0.
    with (out/'gradient_comparisons.jsonl').open('w') as handle:
        for ix,entry in enumerate(entries):
            assert sha(entry['path'])==entry['sha256'];saved=json.loads(Path(entry['path']).read_text());assert saved['plan_sha256']==ph
            payload=saved['payload'];assert payload['status']!='fit_error_requires_review'
            bg,fam,indices,x,y=prepared[entry['fit_input_id']]
            cache=CachedMatchedLikelihood(bg,fam,factors[entry['tree']][indices],x,y)
            theta=np.array(payload['log1p_ratios']);upper=np.log1p(payload['maximum_ratio'])
            analytic=evaluate_gradient(cache,np.expm1(theta));g=analytic['log1p_ratio_gradient']
            np.testing.assert_allclose(analytic['negative_profiled_reml'],payload['negative_profiled_reml'],rtol=1e-9,atol=1e-7)
            max_objective_error=max(max_objective_error,abs(analytic['negative_profiled_reml']-payload['negative_profiled_reml']))
            projected=np.where(theta<=1e-7,np.minimum(g,0),np.where(theta>=upper-1e-7,np.maximum(g,0),g))
            finite={}
            for step in [1e-3,1e-4,1e-5,1e-6]:
                values=[]
                for k in range(3):
                    lo=theta.copy();hi=theta.copy();lo[k]=max(0.,theta[k]-step);hi[k]=min(upper,theta[k]+step)
                    values.append((cache.evaluate(*np.expm1(hi))['negative_profiled_reml']-cache.evaluate(*np.expm1(lo))['negative_profiled_reml'])/(hi[k]-lo[k]))
                finite[str(step)]=values
            np.testing.assert_allclose(finite['0.0001'],payload['finite_difference_gradient'],rtol=1e-4,atol=1e-5)
            below=bool(np.max(abs(projected))<=1e-3)
            counts[f"original_gradient_pass={payload['checks']['projected_gradient_pass']},analytic_below_same_threshold={below}"]+=1
            row=dict(fit_input_id=entry['fit_input_id'],tree=entry['tree'],original_status=payload['status'],
                     analytic_log1p_gradient=g.tolist(),analytic_projected_gradient=projected.tolist(),
                     analytic_below_original_threshold=below,finite_difference_gradients=finite)
            handle.write(json.dumps(row)+'\n')
            if (ix+1)%20==0:print('Diagnosed frozen-prefix gradients',ix+1,'/ 320',flush=True)
    receipt=dict(status='complete_frozen_prefix_gradient_diagnostic',fits=len(entries),comparisons=dict(counts),
                 maximum_objective_error=max_objective_error,snapshot_sha256=sha(snapshot_path),plan_sha256=ph,gradient_proof_sha256=sha(proof_path),
                 script_sha256=sha(__file__),artifacts={'gradient_comparisons.jsonl':sha(out/'gradient_comparisons.jsonl')},
                 scope='All 320 fits in the preexisting frozen prefix; no favorable subset. Analytic score independently validated on dense fixtures, compared with four finite-difference steps and original saved objective/gradient. Original statuses remain unchanged. Not a full-grid reclassification, global-optimum proof or biological inference.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)


if __name__=='__main__':main()
