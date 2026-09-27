"""Assess every exact unique production input for joint zero-reference support."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
import psutil
import subprocess
import time
from matched_joint_covariate_support import assess_support
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
    dependency = json.loads(Path(config['inventory_audit_launch']).read_text())
    while True:
        try:
            process = psutil.Process(dependency['pid'])
            if process.create_time() != dependency['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            assert process.cmdline() == dependency['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',dependency['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state == dict(ActiveState='inactive',Result='success',ExecMainStatus='0'), state
    verify()
    plan = json.loads(Path(config['inventory_plan']).read_text())
    inventory = Path(plan['output'])
    receipt = json.loads((inventory/'receipt.json').read_text())
    audit = json.loads(Path(config['inventory_audit']).read_text())
    assert audit['status'] == 'passed_full_nonlinear_model_input_inventory_readback'
    assert audit['source_receipt_sha256'] == sha(inventory/'receipt.json')
    recipes = [json.loads(line) for line in (inventory/'unique_fit_recipes.jsonl').open()]
    assert len(recipes) == receipt['unique_record_inputs'] == 57616
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
    with (out/'joint_support.jsonl').open('x') as handle:
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
                    result = assess_support(matrix[:,1:])
                    result.update(polynomial_degree=degree,fit_input_id=identifier,partition_index=part['index'],values_sha256=recipe['values_sha256'],ordered_identity_sha256=recipe['ordered_identity_sha256'])
                    handle.write(json.dumps(result,allow_nan=False)+'\n')
                    counts[result['classification']] += 1
                    seen.add(identifier)
                    if len(seen)%500 == 0:
                        handle.flush()
                        print('Joint-support inputs',len(seen),'/ 57616',dict(counts),flush=True)
    assert seen == {r['fit_input_id'] for r in recipes}
    verify()
    result = dict(status='complete_full_nonlinear_joint_support_pending_readback',
        unique_inputs=len(seen),classification_counts=dict(counts),plan_sha256=ch,
        source_inventory_receipt_sha256=sha(inventory/'receipt.json'),scipy_version=scipy.__version__,
        artifacts={'joint_support.jsonl':sha(out/'joint_support.jsonl')},
        scope='Every exact quadratic/cubic candidate input retained with primal support or separating-vector certificate where resolved. Positive diagonal scaling without centering; constant covariates retained. Numerical hull support is not interior overlap, dense local sampling, model adequacy, causal exchangeability or calibrated inference. Same covariates reused by all five trees. Original fits and statuses unchanged; full serialized certificate readback pending.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
