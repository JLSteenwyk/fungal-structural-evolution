"""Reconstruct every exact original fit input once for subsequent simulation."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha, write_json


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    config=json.loads(args.plan.read_text());ch=sha(args.plan)
    production_path=Path(config['production_plan']);plan=json.loads(production_path.read_text())
    def verify():
        assert sha(args.plan)==ch
        for path,h in {**plan['pins'],**config['pins']}.items():assert sha(path)==h,path
    verify();threadpool_limits(1)
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=True)
    lock=(out/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (out/'receipt.json').exists():raise FileExistsError('Cache already completed')
    snapshot=out/'run_plan.json'
    if snapshot.exists():assert snapshot.read_bytes()==args.plan.read_bytes()
    else:snapshot.write_bytes(args.plan.read_bytes())
    inventory=Path(plan['inventory'])
    proof=json.loads(Path(plan['inventory_audit']).read_text())
    assert proof['status']=='passed_full_matched_fit_inventory_readback'
    assert proof['source_receipt_sha256']==sha(inventory/'receipt.json')
    recipes=[json.loads(line) for line in (inventory/'unique_fit_recipes.jsonl').open()]
    assert len(recipes)==len({r['fit_input_id'] for r in recipes})==28808
    assert sum(r['records'] for r in recipes)==130910712
    factors=Path(plan['factors']);factor_files=list(factors.glob('*.npz'));assert len(factor_files)==5
    factor_dir=out/'factors';factor_dir.mkdir(exist_ok=True)
    factor_shapes={}
    for path in factor_files:
        target=factor_dir/path.name
        if target.exists():assert sha(target)==sha(path)
        else:shutil.copy2(path,target)
        with np.load(target,allow_pickle=False) as arrays:factor_shapes[path.stem]=arrays['factor'].shape
    pattern=pd.read_csv(factors/'patterns.tsv',sep='\t').set_index('species_pattern_id').row_index
    nodes=pd.read_csv(plan['nodes'],sep='\t')
    targets=nodes[nodes.role.eq('target')][['node_id','guide','family_component']].rename(columns={'node_id':'target_id'})
    pairs=pd.read_csv(plan['pairs'],sep='\t',usecols=['target_id','background_id','species_pattern_id']).merge(targets,on='target_id',validate='many_to_one')
    pairs['row_identity']=[hashlib.sha256(json.dumps(list(r),separators=(',',':')).encode()).hexdigest() for r in pairs[['target_id','background_id','family_component','species_pattern_id']].itertuples(index=False,name=None)]
    selections=pd.read_csv(plan['selections'],sep='\t',usecols=['target_id','background_id','domain_config_id','policy','scenario_id']).merge(pairs,on=['target_id','background_id'],validate='many_to_one').sort_values('target_id',kind='stable')
    assert len(selections)==2786912
    source=Path(plan['summaries']);numeric=['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference']
    seen=set();manifest=[]
    for part in json.loads((source/'partition_manifest.json').read_text()):
        wanted={(r['guide'],r['policy'],r['scenario_id']):r for r in recipes if r['partition_index']==part['index']}
        if not wanted:continue
        if shutil.disk_usage(out).free<config['resources']['minimum_free_disk_gib']*2**30:raise RuntimeError('Disk reserve reached')
        path=source/part['path'];assert sha(path)==part['sha256']
        records=selections.merge(pd.read_parquet(path,columns=['domain_config_id']+numeric),on='domain_config_id',validate='many_to_one',sort=False)
        for key,frame in records.groupby(['guide','policy','scenario_id'],sort=True):
            if key not in wanted:continue
            recipe=wanted[key];identifier=recipe['fit_input_id'];assert identifier not in seen;seen.add(identifier)
            assert frame.target_id.is_monotonic_increasing and len(frame)==recipe['records']
            matrix=np.ascontiguousarray(frame[numeric].to_numpy(),dtype='<f8');matrix[matrix==0]=0.
            identities=np.asarray(frame.row_identity,dtype='S64')
            assert np.isfinite(matrix).all() and hashlib.sha256(matrix.tobytes()).hexdigest()==recipe['values_sha256']
            assert hashlib.sha256(identities.tobytes()).hexdigest()==recipe['ordered_identity_sha256']
            indices=frame.species_pattern_id.map(pattern);assert indices.notna().all()
            indices=indices.to_numpy(dtype='<i8')
            assert indices.min()>=0 and all(indices.max()<shape[0] for shape in factor_shapes.values())
            active=np.ptp(matrix[:,1:],axis=0)>1e-12
            assert np.all(abs(matrix[:,1:][:,~active])<=1e-12)
            arrays=dict(matrix=matrix,row_identity=identities,
                background=pd.factorize(frame.background_id,sort=True)[0].astype('<i8'),
                family=pd.factorize(frame.family_component,sort=True)[0].astype('<i8'),
                pattern_rows=indices,active_covariates=active,
                covariate_scales=np.std(matrix[:,1:][:,active],axis=0))
            target=out/'inputs'/identifier[:2]/(identifier+'.npz');target.parent.mkdir(parents=True,exist_ok=True)
            if not target.exists():
                temporary=target.with_suffix('.tmp')
                with temporary.open('wb') as f:np.savez_compressed(f,**arrays)
                temporary.replace(target)
            with np.load(target,allow_pickle=False) as saved:
                assert set(saved.files)==set(arrays)
                for name,value in arrays.items():np.testing.assert_array_equal(saved[name],value)
            manifest.append(dict(fit_input_id=identifier,path=str(target.relative_to(out)),sha256=sha(target),records=len(frame),recipe=recipe))
        print('Cached exact unique inputs',len(seen),'/28808',flush=True)
    assert seen=={r['fit_input_id'] for r in recipes}
    verify()
    (out/'manifest.jsonl').write_text(''.join(json.dumps(r,sort_keys=True)+'\n' for r in sorted(manifest,key=lambda x:x['fit_input_id'])))
    write_json(out/'receipt.json',dict(status='all_exact_simulation_inputs_cached_pending_likelihood_replay',
        plan_sha256=ch,unique_inputs=len(seen),record_instances=sum(r['records'] for r in manifest),trees=5,
        artifacts={str(p.relative_to(out)):sha(p) for p in [out/'manifest.jsonl',snapshot,*factor_dir.glob('*.npz')]},
        scope='Exact source numeric/identity hashes, all five factor bindings and full array serialization checked. No row subsampling. Independent likelihood replay from cache, full simulation, calibrated uncertainty and adequacy remain pending.'))


if __name__=='__main__':main()
