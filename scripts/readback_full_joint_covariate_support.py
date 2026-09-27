"""Read back every joint-support certificate and map it to the complete settings grid."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
import subprocess
import time
import psutil
from check_joint_support_certificate import check_certificate
from screen_duplication_domain_alignment_coverage import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--launch', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.plan.read_text())
    ch = sha(args.plan)
    def verify():
        assert sha(args.plan) == ch
        for path, digest in config['pins'].items():
            assert sha(path) == digest, path
    verify()
    launch_hash = sha(args.launch)
    launch = json.loads(args.launch.read_text())
    assert launch['plan_sha256'] == ch
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
    assert sha(args.launch) == launch_hash
    verify()
    produced = Path(config['output'])
    support_hash = sha(produced/'receipt.json')
    support_receipt = json.loads((produced/'receipt.json').read_text())
    assert support_receipt['status'] == 'complete_full_joint_covariate_support_pending_readback'
    assert support_receipt['plan_sha256'] == ch
    assert sha(produced/'joint_support.jsonl') == support_receipt['artifacts']['joint_support.jsonl']
    solutions = {}
    for line in (produced/'joint_support.jsonl').open():
        saved = json.loads(line)
        assert saved['fit_input_id'] not in solutions
        solutions[saved['fit_input_id']] = saved
    assert len(solutions) == support_receipt['unique_inputs'] == 28808
    plan = json.loads(Path(config['production_plan']).read_text())
    inventory = Path(plan['inventory'])
    receipt = json.loads((inventory/'receipt.json').read_text())
    audit = json.loads(Path(plan['inventory_audit']).read_text())
    assert audit['status'] == 'passed_full_matched_fit_inventory_readback'
    assert audit['source_receipt_sha256'] == sha(inventory/'receipt.json')
    recipes = [json.loads(line) for line in (inventory/'unique_fit_recipes.jsonl').open()]
    assert len(recipes) == receipt['unique_record_inputs'] == 28808
    nodes = pd.read_csv(plan['nodes'], sep='\t')
    targets = nodes[nodes.role.eq('target')][['node_id','guide','family_component']].rename(columns={'node_id':'target_id'})
    pairs = pd.read_csv(plan['pairs'],sep='\t',usecols=['target_id','background_id','species_pattern_id']).merge(targets,on='target_id',validate='many_to_one')
    pairs['row_identity'] = [hashlib.sha256(json.dumps(list(r),separators=(',',':')).encode()).hexdigest() for r in pairs[['target_id','background_id','family_component','species_pattern_id']].itertuples(index=False,name=None)]
    selected = pd.read_csv(plan['selections'],sep='\t',usecols=['target_id','background_id','domain_config_id','policy','scenario_id']).merge(pairs,on=['target_id','background_id'],validate='many_to_one').sort_values('target_id',kind='stable')
    assert len(selected) == 2786912
    numeric = ['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference']
    source = Path(plan['summaries'])
    out = args.output
    out.mkdir(parents=True,exist_ok=False)
    counts, seen = Counter(), set()
    with (out/'joint_support.jsonl').open('x') as handle:
        for part in json.loads((source/'partition_manifest.json').read_text()):
            wanted = {(r['guide'],r['policy'],r['scenario_id']): r for r in recipes if r['partition_index']==part['index']}
            if not wanted:
                continue
            path = source/part['path']
            assert sha(path) == part['sha256']
            records = selected.merge(pd.read_parquet(path,columns=['domain_config_id']+numeric),on='domain_config_id',validate='many_to_one',sort=False)
            for key, frame in records.groupby(['guide','policy','scenario_id'],sort=True):
                if key not in wanted:
                    continue
                recipe = wanted[key]
                identifier = recipe['fit_input_id']
                assert identifier not in seen and frame.target_id.is_monotonic_increasing
                matrix = np.ascontiguousarray(frame[numeric].to_numpy(),dtype='<f8')
                matrix[matrix==0] = 0.
                assert len(frame) == recipe['records']
                assert hashlib.sha256(matrix.tobytes()).hexdigest() == recipe['values_sha256']
                assert hashlib.sha256(np.asarray(frame.row_identity,dtype='S64').tobytes()).hexdigest() == recipe['ordered_identity_sha256']
                result = solutions[identifier]
                assert result['partition_index'] == part['index']
                assert result['values_sha256'] == recipe['values_sha256']
                assert result['ordered_identity_sha256'] == recipe['ordered_identity_sha256']
                check_certificate(matrix[:,1:], result)
                result.update(fit_input_id=identifier,partition_index=part['index'],values_sha256=recipe['values_sha256'],ordered_identity_sha256=recipe['ordered_identity_sha256'])
                handle.write(json.dumps(result,allow_nan=False)+'\n')
                counts[result['classification']] += 1
                seen.add(identifier)
                if len(seen)%500 == 0:
                    handle.flush()
                    print('Joint-support inputs',len(seen),'/ 28808',dict(counts),flush=True)
    assert seen == {r['fit_input_id'] for r in recipes}
    assert seen == set(solutions)
    assert dict(counts) == support_receipt['classification_counts']
    assert sha(produced/'receipt.json') == support_hash
    assert sha(produced/'joint_support.jsonl') == support_receipt['artifacts']['joint_support.jsonl']
    map_path = inventory/'full_setting_fit_map.tsv'
    assert sha(map_path) == receipt['artifacts'][map_path.name]
    mapping = pd.read_csv(map_path,sep='\t')
    assert len(mapping) == 82944 and set(mapping.fit_input_id) == seen
    mapping['joint_support_classification'] = mapping.fit_input_id.map({k:v['classification'] for k,v in solutions.items()})
    mapping.to_csv(out/'full_setting_support.tsv',sep='\t',index=False)
    exported = pd.read_csv(out/'full_setting_support.tsv',sep='\t')
    pd.testing.assert_frame_equal(exported.drop(columns='joint_support_classification'),pd.read_csv(map_path,sep='\t'))
    assert exported.joint_support_classification.tolist() == [solutions[k]['classification'] for k in exported.fit_input_id]
    verify()
    result = dict(status='passed_full_joint_covariate_support_certificate_readback',
        source_support_receipt_sha256=support_hash, full_settings=len(mapping), setting_classification_counts=mapping.joint_support_classification.value_counts().to_dict(),
        script_sha256=sha(__file__), checker_sha256=sha(Path(__file__).with_name('check_joint_support_certificate.py')), unique_inputs=len(seen),classification_counts=dict(counts),plan_sha256=ch,
        source_inventory_receipt_sha256=sha(inventory/'receipt.json'),scipy_version=scipy.__version__,
        artifacts={name:sha(out/name) for name in ['joint_support.jsonl','full_setting_support.tsv']},
        scope='Every exact production input retained with primal support or separating-vector certificate where resolved. Positive diagonal scaling without centering; constant covariates retained. Numerical hull support is not interior overlap, dense local sampling, model adequacy, causal exchangeability or calibrated inference. Same covariates reused by all five trees. Original fits and statuses unchanged; all saved certificates checked without an optimizer, and all original settings mapped/read back. Does not re-solve failed optimization cases.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
