"""Try explicit nonnegative certificates for every unresolved joint-support input."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
import psutil
from check_joint_support_certificate import check_certificate
from project_joint_support_weights import project_weights
from screen_duplication_domain_alignment_coverage import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();launch=json.loads(Path(plan['dependency_launch']).read_text())
    while True:
        try:
            process=psutil.Process(launch['pid'])
            if process.create_time()!=launch['created'] or process.status()==psutil.STATUS_ZOMBIE:break
            assert process.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True)
    assert dict(l.split('=',1) for l in state.splitlines())==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    verify();root=Path(plan['readback']);rp=root/'receipt.json';rh=sha(rp);r=json.loads(rp.read_text())
    assert r['status']=='passed_full_joint_covariate_support_certificate_readback'
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    saved={}
    for line in (root/'joint_support.jsonl').open():
        row=json.loads(line);assert row['fit_input_id'] not in saved;saved[row['fit_input_id']]=row
    assert len(saved)==28808
    wanted={k:v for k,v in saved.items() if v['classification'].startswith('unresolved')}
    production=json.loads(Path(plan['production_plan']).read_text());inventory=Path(production['inventory'])
    recipes=[json.loads(l) for l in (inventory/'unique_fit_recipes.jsonl').open()]
    recipes=[x for x in recipes if x['fit_input_id'] in wanted]
    assert {x['fit_input_id'] for x in recipes}==set(wanted)
    nodes=pd.read_csv(production['nodes'],sep='\t');targets=nodes[nodes.role.eq('target')][['node_id','guide','family_component']].rename(columns={'node_id':'target_id'})
    pairs=pd.read_csv(production['pairs'],sep='\t',usecols=['target_id','background_id','species_pattern_id']).merge(targets,on='target_id',validate='many_to_one')
    pairs['row_identity']=[hashlib.sha256(json.dumps(list(x),separators=(',',':')).encode()).hexdigest() for x in pairs[['target_id','background_id','family_component','species_pattern_id']].itertuples(index=False,name=None)]
    selected=pd.read_csv(production['selections'],sep='\t',usecols=['target_id','background_id','domain_config_id','policy','scenario_id']).merge(pairs,on=['target_id','background_id'],validate='many_to_one').sort_values('target_id',kind='stable')
    assert len(selected)==2786912
    numeric=['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference']
    source=Path(production['summaries']);results={}
    for part in json.loads((source/'partition_manifest.json').read_text()):
        group={(x['guide'],x['policy'],x['scenario_id']):x for x in recipes if x['partition_index']==part['index']}
        if not group:continue
        path=source/part['path'];assert sha(path)==part['sha256']
        records=selected.merge(pd.read_parquet(path,columns=['domain_config_id']+numeric),on='domain_config_id',validate='many_to_one',sort=False)
        for key,frame in records.groupby(['guide','policy','scenario_id'],sort=True):
            if key not in group:continue
            recipe=group[key];identifier=recipe['fit_input_id'];assert identifier not in results
            matrix=np.ascontiguousarray(frame[numeric].to_numpy(),dtype='<f8');matrix[matrix==0]=0.
            assert len(frame)==recipe['records'] and frame.target_id.is_monotonic_increasing
            assert hashlib.sha256(matrix.tobytes()).hexdigest()==recipe['values_sha256']
            assert hashlib.sha256(np.asarray(frame.row_identity,dtype='S64').tobytes()).hexdigest()==recipe['ordered_identity_sha256']
            original=wanted[identifier];check_certificate(matrix[:,1:],original)
            result=project_weights(matrix[:,1:],original)
            result.update(fit_input_id=identifier,original_classification=original['classification'],values_sha256=recipe['values_sha256'],ordered_identity_sha256=recipe['ordered_identity_sha256'])
            results[identifier]=result
        print('Projected certificates',len(results),'/',len(wanted),flush=True)
    assert set(results)==set(wanted)
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    path=out/'projected_certificates.jsonl';path.write_text(''.join(json.dumps(v,allow_nan=False)+'\n' for v in results.values()))
    assert {v['fit_input_id']:v for v in map(json.loads,path.read_text().splitlines())}==results
    mapping=pd.read_csv(root/'full_setting_support.tsv',sep='\t')
    mapping['projected_support_classification']=mapping.fit_input_id.map({k:v['classification'] for k,v in results.items()}).fillna('original_certificate_retained')
    mapping.to_csv(out/'full_setting_projection.tsv',sep='\t',index=False)
    pd.testing.assert_frame_equal(pd.read_csv(out/'full_setting_projection.tsv',sep='\t'),mapping)
    verify();assert sha(rp)==rh
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    receipt=dict(status='complete_unresolved_joint_support_projection',plan_sha256=ph,source_readback_receipt_sha256=rh,original_inputs=len(saved),unresolved_inputs=len(wanted),classification_counts=dict(Counter(v['classification'] for v in results.values())),full_settings=len(mapping),setting_classification_counts=mapping.projected_support_classification.value_counts().to_dict(),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All unresolved inputs retained. Nonnegative projected weights checked against original fingerprinted matrices using scalar and vector barycenters. Numerical zero support only; no LP optimum, interior overlap, causal adequacy or change to original classifications/fit flags.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':main()
