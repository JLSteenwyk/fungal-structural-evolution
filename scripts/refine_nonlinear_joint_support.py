"""Attempt explicit nonnegative certificates for all unresolved nonlinear inputs, preserving original classifications."""
import argparse
from collections import Counter
import hashlib
import math
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
import subprocess
import time
import psutil
from check_joint_support_certificate import check_certificate
from project_joint_support_weights import project_weights
from screen_duplication_domain_alignment_coverage import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.plan.read_text())
    ch = sha(args.plan)
    def verify():
        assert sha(args.plan) == ch
        for path, digest in config['pins'].items():
            assert sha(path) == digest, path
    verify()
    launch_hash = sha(config['dependency_launch'])
    launch = json.loads(Path(config['dependency_launch']).read_text())
    assert launch['plan_sha256'] == sha(config['production_plan'])
    while True:
        try:
            process = psutil.Process(launch['pid'])
            if process.create_time() != launch['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            assert process.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    raw = subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True)
    assert dict(line.split('=',1) for line in raw.splitlines()) == dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    assert sha(config['dependency_launch']) == launch_hash
    verify()
    produced = Path(config['readback'])
    support_hash = sha(produced/'receipt.json')
    support_receipt = json.loads((produced/'receipt.json').read_text())
    assert support_receipt['status'] == 'passed_full_nonlinear_joint_support_certificate_readback'
    assert support_receipt['plan_sha256'] == sha(config['production_plan'])
    assert sha(produced/'joint_support.jsonl') == support_receipt['artifacts']['joint_support.jsonl']
    for name,digest in support_receipt['artifacts'].items():assert sha(produced/name)==digest
    solutions = {}
    for line in (produced/'joint_support.jsonl').open():
        saved = json.loads(line)
        assert saved['fit_input_id'] not in solutions
        solutions[saved['fit_input_id']] = saved
    assert len(solutions) == support_receipt['unique_inputs'] == 57616
    production = json.loads(Path(config['production_plan']).read_text())
    plan = json.loads(Path(production['inventory_plan']).read_text())
    inventory = Path(plan['output'])
    receipt = json.loads((inventory/'receipt.json').read_text())
    audit = json.loads(Path(production['inventory_audit']).read_text())
    assert audit['status'] == 'passed_full_nonlinear_model_input_inventory_readback'
    assert audit['source_receipt_sha256'] == sha(inventory/'receipt.json')
    recipes = [json.loads(line) for line in (inventory/'unique_fit_recipes.jsonl').open()]
    assert len(recipes) == receipt['unique_record_inputs'] == 57616
    unresolved = {k:v for k,v in solutions.items() if v['classification'].startswith('unresolved')}
    recipes = [r for r in recipes if r['fit_input_id'] in unresolved]
    assert {r['fit_input_id'] for r in recipes} == set(unresolved)
    nodes = pd.read_csv(plan['nodes'], sep='\t')
    targets = nodes[nodes.role.eq('target')][['node_id','guide','family_component']].rename(columns={'node_id':'target_id'})
    pairs = pd.read_csv(plan['pairs'],sep='\t',usecols=['target_id','background_id','species_pattern_id']).merge(targets,on='target_id',validate='many_to_one')
    pairs['row_identity'] = [hashlib.sha256(json.dumps(list(r),separators=(',',':')).encode()).hexdigest() for r in pairs[['target_id','background_id','family_component','species_pattern_id']].itertuples(index=False,name=None)]
    selected = pd.read_csv(plan['selections'],sep='\t',usecols=['target_id','background_id','domain_config_id','policy','scenario_id']).merge(pairs,on=['target_id','background_id'],validate='many_to_one').sort_values('target_id',kind='stable')
    assert len(selected) == 2786912
    all_numeric = ['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference','identity_power_2_difference','identity_power_3_difference']
    nonlinear = Path(plan['nonlinear'])
    nr = json.loads((nonlinear/'receipt.json').read_text())
    for name,digest in nr['artifacts'].items():assert sha(nonlinear/name)==digest
    source = Path(plan['summaries'])
    out = Path(config['output'])
    out.mkdir(parents=True,exist_ok=False)
    counts, seen = Counter(), set()
    projections = {}
    with (out/'projected_certificates.jsonl').open('x') as handle:
        for part in json.loads((source/'partition_manifest.json').read_text()):
            wanted = {(r['guide'],r['policy'],r['scenario_id'],r['polynomial_degree']): r for r in recipes if r['partition_index']==part['index']}
            if not wanted:
                continue
            path = source/part['path']
            assert sha(path) == part['sha256']
            values = pd.read_parquet(path,columns=['domain_config_id']+all_numeric[:5])
            extra = pd.read_parquet(nonlinear/f"{part['index']:03d}.parquet")
            assert set(values.domain_config_id)==set(extra.domain_config_id)
            values = values.merge(extra[['domain_config_id']+all_numeric[5:]],on='domain_config_id',validate='one_to_one')
            records = selected.merge(values,on='domain_config_id',validate='many_to_one',sort=False)
            for key, frame in records.groupby(['guide','policy','scenario_id'],sort=True):
                for degree in [2,3]:
                    numeric = all_numeric[:degree+4]
                    if (*key,degree) not in wanted:
                        continue
                    recipe = wanted[(*key,degree)]
                    identifier = recipe['fit_input_id']
                    assert identifier not in seen and frame.target_id.is_monotonic_increasing
                    matrix = np.ascontiguousarray(frame[numeric].to_numpy(),dtype='<f8')
                    matrix[matrix==0] = 0.
                    assert len(frame) == recipe['records']
                    assert hashlib.sha256(matrix.tobytes()).hexdigest() == recipe['values_sha256']
                    assert hashlib.sha256(np.asarray(frame.row_identity,dtype='S64').tobytes()).hexdigest() == recipe['ordered_identity_sha256']
                    result = solutions[identifier]
                    assert result['polynomial_degree'] == degree
                    assert result['partition_index'] == part['index']
                    assert result['values_sha256'] == recipe['values_sha256']
                    assert result['ordered_identity_sha256'] == recipe['ordered_identity_sha256']
                    check_certificate(matrix[:,1:], result)
                    original = result
                    result = project_weights(matrix[:,1:], original)
                    result['original_classification'] = original['classification']
                    result.update(polynomial_degree=degree,fit_input_id=identifier,partition_index=part['index'],values_sha256=recipe['values_sha256'],ordered_identity_sha256=recipe['ordered_identity_sha256'])
                    serialized = json.dumps(result,allow_nan=False)
                    restored = json.loads(serialized)
                    assert restored == result
                    if 'support_weights' in restored:
                        weights = restored['support_weights']; indices = restored['support_indices']
                        assert len(weights)==len(indices) and len(set(indices))==len(indices)
                        assert all(math.isfinite(w) and w>=0 for w in weights)
                        assert all(type(i) is int and 0<=i<len(matrix) for i in indices)
                        assert abs(math.fsum(weights)-1)<=1e-12
                        x=matrix[:,1:];scales=np.array([max(abs(float(v)) for v in x[:,j]) or 1. for j in range(x.shape[1])])
                        np.testing.assert_array_equal(scales,restored['scales'])
                        bary=[math.fsum(w*float(x[i,j]/scales[j]) for i,w in zip(indices,weights)) for j in range(x.shape[1])]
                        np.testing.assert_allclose(bary,restored['barycenter'],rtol=1e-10,atol=1e-12)
                        distance=max(map(abs,bary));assert math.isclose(distance,restored['primal_distance'],rel_tol=1e-10,abs_tol=1e-12)
                        expected='supported_by_projected_nonnegative_weights' if distance<=1e-8 else 'unresolved_after_projection'
                        assert restored['classification']==expected
                    else:
                        assert restored['classification']=='unresolved_no_positive_weights'
                        assert not any(w>0 for w in original.get('support_weights',[]))
                    handle.write(serialized+'\n')
                    projections[identifier]=restored
                    counts[result['classification']] += 1
                    seen.add(identifier)
                    if len(seen)%500 == 0:
                        handle.flush()
                        print('Projected unresolved inputs',len(seen),'/',len(unresolved),dict(counts),flush=True)
    assert seen == {r['fit_input_id'] for r in recipes}
    assert seen == set(unresolved)
    assert {v['fit_input_id']:v for v in map(json.loads,(out/'projected_certificates.jsonl').read_text().splitlines())} == projections
    assert sha(produced/'receipt.json') == support_hash
    assert sha(produced/'joint_support.jsonl') == support_receipt['artifacts']['joint_support.jsonl']
    map_path = inventory/'full_setting_fit_map.tsv'
    assert sha(map_path) == receipt['artifacts'][map_path.name]
    mapping = pd.read_csv(map_path,sep='\t')
    assert len(mapping) == 165888 and set(mapping.fit_input_id) == set(solutions)
    mapping['joint_support_classification'] = mapping.fit_input_id.map({k:v['classification'] for k,v in solutions.items()})
    mapping['projected_support_classification']=mapping.fit_input_id.map({k:v['classification'] for k,v in projections.items()}).fillna('original_certificate_retained')
    mapping.to_csv(out/'full_setting_support.tsv',sep='\t',index=False)
    exported = pd.read_csv(out/'full_setting_support.tsv',sep='\t')
    pd.testing.assert_frame_equal(exported.drop(columns=['joint_support_classification','projected_support_classification']),pd.read_csv(map_path,sep='\t'))
    assert exported.joint_support_classification.tolist() == [solutions[k]['classification'] for k in exported.fit_input_id]
    assert exported.projected_support_classification.tolist()==[projections[k]['classification'] if k in projections else 'original_certificate_retained' for k in exported.fit_input_id]
    verify()
    result = dict(status='complete_nonlinear_unresolved_joint_support_projection',
        dependency_terminal_state=dict(line.split('=',1) for line in raw.splitlines()),
        original_inputs=len(solutions),unresolved_inputs=len(unresolved),projected_setting_classification_counts=mapping.projected_support_classification.value_counts().to_dict(),
        source_support_receipt_sha256=support_hash, full_settings=len(mapping), setting_classification_counts=mapping.joint_support_classification.value_counts().to_dict(),
        script_sha256=sha(__file__), checker_sha256=sha(Path(__file__).with_name('check_joint_support_certificate.py')), unique_inputs=len(seen),classification_counts=dict(counts),plan_sha256=ch,
        source_inventory_receipt_sha256=sha(inventory/'receipt.json'),scipy_version=scipy.__version__,
        artifacts={name:sha(out/name) for name in ['projected_certificates.jsonl','full_setting_support.tsv']},
        scope='All unresolved quadratic/cubic inputs retained, with separate nonnegative projected weight certificates checked against original fingerprinted matrices and serialized outputs. Original classifications remain unchanged in the full setting map. No LP optimum, interior overlap, model adequacy or calibrated biological inference claimed.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
