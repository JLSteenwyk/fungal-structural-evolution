"""Replay every original fit from the complete serialized simulation-input cache."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import psutil
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha, write_json
from matched_mixed_covariance import MatchedCovariance, profiled_reml


def replay(arrays, old, factor):
    matrix=arrays['matrix'];n=len(matrix)
    assert matrix.shape==(n,5) and np.isfinite(matrix).all() and n==old['records']
    active=np.ptp(matrix[:,1:],axis=0)>1e-12
    np.testing.assert_array_equal(active,arrays['active_covariates'])
    np.testing.assert_array_equal(active,old['active_covariates'])
    assert np.all(abs(matrix[:,1:][:,~active])<=1e-12)
    scales=np.std(matrix[:,1:][:,active],axis=0)
    np.testing.assert_array_equal(scales,arrays['covariate_scales'])
    np.testing.assert_allclose(scales,old['covariate_scales'],rtol=1e-13,atol=1e-14)
    indices=arrays['pattern_rows'];assert indices.shape==(n,) and indices.dtype.kind in 'iu'
    assert indices.min()>=0 and indices.max()<factor.shape[0]
    x=np.column_stack([np.ones(n),matrix[:,1:][:,active]/scales])
    result=profiled_reml(MatchedCovariance(arrays['background'],arrays['family'],factor[indices],1.,*old['ratios']),x,matrix[:,0])
    for key in ['negative_profiled_reml','beta','profiled_scale','residual_quadratic','conditional_beta_covariance']:
        np.testing.assert_allclose(result[key],old[key],rtol=1e-7,atol=1e-8,err_msg=key)
    conversion=np.r_[1.,1/scales]
    np.testing.assert_allclose(result['beta']*conversion,old['raw_unit_beta'],rtol=1e-7,atol=1e-8)
    np.testing.assert_allclose(result['conditional_beta_covariance']*conversion[:,None]*conversion[None,:],old['raw_unit_conditional_beta_covariance'],rtol=1e-7,atol=1e-8)
    return dict(objective_absolute_error=abs(result['negative_profiled_reml']-old['negative_profiled_reml']),
                beta_maximum_absolute_error=float(np.max(abs(result['beta']-old['beta']))))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    config=json.loads(args.plan.read_text());ch=sha(args.plan);threadpool_limits(1)
    def verify():
        assert sha(args.plan)==ch
        for name,h in config['pins'].items():assert sha(name)==h,name
    verify();launch=json.loads(Path(config['producer_launch']).read_text())
    while True:
        state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','MainPID','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:break
        assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid']
        try:
            p=psutil.Process(launch['pid']);assert p.create_time()==launch['created'] and p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:time.sleep(1);continue
        time.sleep(30)
    assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0'
    verify();cache=Path(config['cache']);source=Path(config['original'])
    cr=json.loads((cache/'receipt.json').read_text());sr=json.loads((source/'receipt.json').read_text())
    assert cr['status']=='all_exact_simulation_inputs_cached_pending_likelihood_replay'
    assert cr['plan_sha256']==sha(config['producer_plan']) and cr['unique_inputs']==28808
    for root,r in [(cache,cr),(source,sr)]:
        for name,h in r['artifacts'].items():assert sha(root/name)==h,name
    original={}
    for line in (source/'fit_manifest.jsonl').open():
        r=json.loads(line);key=r['fit_input_id'],r['tree'];assert key not in original;original[key]=r
    assert len(original)==144040
    factors={p.stem:np.load(p,allow_pickle=False)['factor'] for p in (cache/'factors').glob('*.npz')};assert len(factors)==5
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=True)
    lock=(out/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    manifest=[];seen=set();fits_seen=set()
    for line in (cache/'manifest.jsonl').open():
        entry=json.loads(line);identifier=entry['fit_input_id'];assert identifier not in seen;seen.add(identifier)
        path=cache/entry['path'];assert path.resolve().is_relative_to(cache.resolve()) and sha(path)==entry['sha256']
        bindings={tree:original[(identifier,tree)]['sha256'] for tree in factors}
        target=out/(identifier+'.json')
        for tree in factors:
            r=original[(identifier,tree)];assert sha(r['path'])==r['sha256'];fits_seen.add((identifier,tree))
        if target.exists():
            result=json.loads(target.read_text());assert result['plan_sha256']==ch and result['cache_sha256']==entry['sha256'] and result['original_fit_hashes']==bindings
            assert set(result['replays'])==set(factors)
        else:
            with np.load(path,allow_pickle=False) as handle:arrays={name:handle[name] for name in handle.files}
            assert hashlib.sha256(arrays['matrix'].tobytes()).hexdigest()==entry['recipe']['values_sha256']
            assert hashlib.sha256(arrays['row_identity'].tobytes()).hexdigest()==entry['recipe']['ordered_identity_sha256']
            checks={}
            for tree,factor in factors.items():
                old=json.loads(Path(original[(identifier,tree)]['path']).read_text())
                assert old['fit_input_id']==identifier and old['tree']==tree and old['plan_sha256']==sr['plan_sha256']
                checks[tree]=replay(arrays,old['payload'],factor)
            result=dict(plan_sha256=ch,cache_sha256=entry['sha256'],original_fit_hashes=bindings,replays=checks)
            write_json(target,result)
        manifest.append(dict(fit_input_id=identifier,path=target.name,sha256=sha(target)))
        if len(seen)%100==0:print('Replayed cached inputs',len(seen),'/28808',flush=True)
    assert len(seen)==28808 and fits_seen==set(original)
    verify();(out/'manifest.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in manifest))
    write_json(out/'receipt.json',dict(status='all_144040_original_fits_replayed_from_simulation_cache',plan_sha256=ch,
        inputs=len(seen),fits=len(fits_seen),source_cache_receipt_sha256=sha(cache/'receipt.json'),
        original_receipt_sha256=sha(source/'receipt.json'),artifacts={'manifest.jsonl':sha(out/'manifest.jsonl')},
        scope='All original fixed-parameter likelihoods, coefficients, scales, quadratics and conditional covariance matrices checked from serialized inputs. Same validated direct evaluator as original fit readback; not independent model adequacy or uncertainty calibration.'))


if __name__=='__main__':main()
