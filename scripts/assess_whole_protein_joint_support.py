#!/usr/bin/env python3
"""Assess zero-reference support for every audited whole-protein input."""
import hashlib
import json
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from assess_whole_protein_model_designs_v2 import matrix, KEYS
from inventory_whole_protein_model_inputs import signature
from matched_joint_covariate_support import assess_support
from screen_duplication_alignment_reuse import sha


def main():
    pp=Path('metadata/whole_protein_joint_support_plan_20260929.json')
    plan=json.loads(pp.read_text());bindings={str(pp):sha(pp),**plan['pins']}
    dependency=json.loads(Path(plan['audit_launch']).read_text())
    while psutil.pid_exists(dependency['pid']):
        try:
            p=psutil.Process(dependency['pid'])
            if abs(p.create_time()-dependency['created'])>.01 or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dependency['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',dependency['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    source_plan=json.loads(Path(plan['inventory_plan']).read_text());root=Path(source_plan['output'])
    receipt=json.loads((root/'receipt.json').read_text());proof=json.loads(Path(plan['audit_proof']).read_text())
    assert proof['status']=='passed_full_whole_protein_input_inventory_readback'
    assert proof['source_receipt_sha256']==sha(root/'receipt.json')
    bindings.update(receipt['source_hashes'])
    bindings.update({str(root/n):h for n,h in receipt['artifacts'].items()})
    bindings[str(root/'receipt.json')]=sha(root/'receipt.json');bindings[plan['audit_proof']]=sha(plan['audit_proof'])
    def verify():
        for path,digest in bindings.items():assert sha(path)==digest,path
    verify()
    recipes=defaultdict(list);expected=set()
    for line in (root/'unique_input_recipes.jsonl').open():
        row=json.loads(line);assert row['fit_input_id'] not in expected
        expected.add(row['fit_input_id']);recipes[row['partition_index']].append(row)
    assert len(expected)==receipt['unique_inputs']
    covariance=Path(source_plan['covariance_root'])
    cov=pd.read_csv(covariance/'pair_covariance_index.tsv',sep='\t',usecols=['target_id','background_id','row_identity'])
    covariance_sha=sha(covariance/'receipt.json')
    selected=pd.read_csv(source_plan['selections'],sep='\t',usecols=['target_id','background_id','policy','scenario_id'])
    assert len(selected)==2786912
    records_root=Path(source_plan['records_root'])
    out=Path(plan['output']);out.mkdir(exist_ok=False)
    seen=set();geometries={};counts=Counter()
    with (out/'input_support_map.jsonl').open('w') as mapping,(out/'joint_support.jsonl').open('w') as certificates:
        for part in json.loads((records_root/'partition_manifest.json').read_text()):
            wanted=recipes.get(part['index'],[])
            if not wanted:continue
            assert sha(records_root/part['path'])==part['sha256']
            values=pd.read_parquet(records_root/part['path'])
            records=selected.merge(values,on=['target_id','background_id'],validate='many_to_one').merge(cov,on=['target_id','background_id'],validate='many_to_one').sort_values('target_id',kind='stable')
            groups=records.groupby(KEYS).indices
            for recipe in wanted:
                key=tuple(recipe[k] for k in KEYS)
                frame=records.iloc[groups[key]]
                chosen,x,names=matrix(frame,recipe['variant'])
                spec=recipe['specification'];columns=spec['columns']
                assert columns[:2]==[recipe['outcome'],'intercept']
                x=x[:,[names.index(c) for c in columns[2:]]]
                numeric=np.column_stack([chosen[recipe['outcome']].to_numpy(),np.ones(len(chosen)),x])
                identifier,_,reconstructed=signature(chosen,numeric,columns,recipe['outcome'],covariance_sha)
                assert identifier==recipe['fit_input_id'] and reconstructed==spec and identifier not in seen
                seen.add(identifier)
                raw=np.ascontiguousarray(x,dtype='<f8');raw[raw==0]=0.
                geometry_spec=dict(records=len(chosen),columns=columns[2:],values_sha256=hashlib.sha256(raw.tobytes()).hexdigest(),ordered_identity_sha256=spec['ordered_identity_sha256'])
                geometry_id=hashlib.sha256(json.dumps(geometry_spec,sort_keys=True,separators=(',',':')).encode()).hexdigest()
                if geometry_id not in geometries:
                    result=assess_support(raw)
                    geometries[geometry_id]=dict(specification=geometry_spec,classification=result['classification'])
                    certificates.write(json.dumps(dict(geometry_id=geometry_id,representative_input=identifier,specification=geometry_spec,result=result),allow_nan=False)+'\n')
                    counts[result['classification']]+=1
                else:assert geometries[geometry_id]['specification']==geometry_spec
                mapping.write(json.dumps(dict(fit_input_id=identifier,geometry_id=geometry_id,classification=geometries[geometry_id]['classification']))+'\n')
                if len(seen)%500==0:
                    mapping.flush();certificates.flush()
                    print('Whole-protein joint support',len(seen),'/',len(expected),'unique geometries',len(geometries),flush=True)
    assert seen==expected
    verify()
    result=dict(status='complete_whole_protein_joint_support_pending_readback',unique_inputs=len(seen),unique_geometries=len(geometries),geometry_classification_counts=dict(counts),source_hashes=bindings,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All audited unique inputs mapped to exact covariate geometry; outcome and tree do not alter hull support. Ordered row identities retained for sparse certificate indices. Positive diagonal scaling, no centering. All unresolved cases preserved. Full certificate readback pending. Hull inclusion does not establish interior overlap, dense support, model adequacy or causal exchangeability.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)


if __name__=='__main__':main()
