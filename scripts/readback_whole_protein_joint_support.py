#!/usr/bin/env python3
"""Reconstruct all geometries and verify support certificates without optimization."""
import hashlib
import json
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from readback_whole_protein_input_inventory import reconstruct, KEYS
from check_joint_support_certificate import check_certificate
from screen_duplication_alignment_reuse import sha


def main():
    pp=Path('metadata/whole_protein_joint_support_readback_plan_20260929.json')
    plan=json.loads(pp.read_text());bindings={str(pp):sha(pp),**plan['pins']}
    launch=json.loads(Path(plan['launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            p=psutil.Process(launch['pid'])
            if abs(p.create_time()-launch['created'])>.01 or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    producer=json.loads(Path(launch['plan']).read_text());root=Path(producer['output'])
    receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='complete_whole_protein_joint_support_pending_readback'
    bindings.update(receipt['source_hashes']);bindings[str(root/'receipt.json')]=sha(root/'receipt.json')
    bindings.update({str(root/n):h for n,h in receipt['artifacts'].items()})
    def verify():
        for path,digest in bindings.items():assert sha(path)==digest,path
    verify()
    source_plan=json.loads(Path(producer['inventory_plan']).read_text());inventory=Path(source_plan['output'])
    recipes=defaultdict(list);expected=set()
    for line in (inventory/'unique_input_recipes.jsonl').open():
        row=json.loads(line);assert row['fit_input_id'] not in expected
        expected.add(row['fit_input_id']);recipes[row['partition_index']].append(row)
    mapping={}
    for line in (root/'input_support_map.jsonl').open():
        row=json.loads(line);assert row['fit_input_id'] not in mapping
        mapping[row['fit_input_id']]=row
    assert set(mapping)==expected and len(expected)==receipt['unique_inputs']
    certificates={}
    for line in (root/'joint_support.jsonl').open():
        row=json.loads(line);assert row['geometry_id'] not in certificates
        certificates[row['geometry_id']]=row
        assert row['representative_input'] in expected
        assert mapping[row['representative_input']]['geometry_id']==row['geometry_id']
    covariance=Path(source_plan['covariance_root']);cov=pd.read_csv(covariance/'pair_covariance_index.tsv',sep='\t')
    identities={}
    for row in cov.to_dict('records'):
        values=[row[k] for k in ['target_id','background_id','family_component','species_pattern_id']]
        digest=hashlib.sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest()
        assert digest==row['row_identity']
        key=row['target_id'],row['background_id'];assert key not in identities
        identities[key]=digest
    assert len(identities)==52675
    selected=pd.read_csv(source_plan['selections'],sep='\t',usecols=['target_id','background_id','policy','scenario_id'])
    assert len(selected)==2786912
    records_root=Path(source_plan['records_root']);seen=set();checked=set();counts=Counter()
    for part in json.loads((records_root/'partition_manifest.json').read_text()):
        wanted=recipes.get(part['index'],[])
        if not wanted:continue
        assert sha(records_root/part['path'])==part['sha256']
        values=pd.read_parquet(records_root/part['path'])
        records=selected.merge(values,on=['target_id','background_id'],validate='many_to_one').sort_values('target_id',kind='stable')
        groups=records.groupby(KEYS).indices
        for recipe in wanted:
            identifier=recipe['fit_input_id'];assert identifier not in seen;seen.add(identifier)
            frame=records.iloc[groups[tuple(recipe[k] for k in KEYS)]]
            chosen,columns=reconstruct(frame,recipe['variant']);spec=recipe['specification']
            x=np.column_stack([columns[n] for n in spec['columns'][2:]]).astype('<f8');x[x==0]=0.
            ordered=b''.join(identities[(t,b)].encode('ascii') for t,b in chosen[['target_id','background_id']].itertuples(index=False,name=None))
            identity_sha=hashlib.sha256(ordered).hexdigest()
            assert identity_sha==spec['ordered_identity_sha256'] and len(chosen)==spec['records']
            numeric=np.column_stack([chosen[recipe['outcome']].to_numpy(),np.ones(len(chosen)),x]).astype('<f8');numeric[numeric==0]=0.
            assert hashlib.sha256(numeric.tobytes()).hexdigest()==spec['values_sha256']
            geometry_spec=dict(records=len(chosen),columns=spec['columns'][2:],values_sha256=hashlib.sha256(x.tobytes()).hexdigest(),ordered_identity_sha256=identity_sha)
            geometry_id=hashlib.sha256(json.dumps(geometry_spec,sort_keys=True,separators=(',',':')).encode()).hexdigest()
            m=mapping[identifier];assert m['geometry_id']==geometry_id
            cert=certificates[geometry_id];assert cert['specification']==geometry_spec
            assert m['classification']==cert['result']['classification']
            if geometry_id not in checked:
                check_certificate(x,cert['result']);checked.add(geometry_id);counts[m['classification']]+=1
        print('Checked joint-support source partition',part['index']+1,'/96',flush=True)
    assert seen==expected and checked==set(certificates) and len(checked)==receipt['unique_geometries']
    assert dict(counts)==receipt['geometry_classification_counts']
    setting_counts=Counter();settings=0
    for line in (inventory/'setting_input_map.jsonl').open():
        row=json.loads(line)
        if row['fit_input_id'] is None:classification='nonestimable_design'
        else:classification=mapping[row['fit_input_id']]['classification']
        setting_counts[classification]+=1;settings+=1
    assert settings==414720
    verify()
    result=dict(status='passed_full_whole_protein_joint_support_readback',unique_inputs=len(seen),unique_geometries=len(checked),geometry_classification_counts=dict(counts),setting_outcome_classification_counts=dict(setting_counts),settings=settings,source_receipt_sha256=sha(root/'receipt.json'),checker_sha256=sha(__file__),producer_terminal_state=state,scope='All covariate geometries and ordered identities rebuilt; every unique certificate checked without optimizer; all setting/outcome links counted. Unresolved statuses retained. No interior-overlap, model adequacy or biological inference claim.')
    with Path(plan['proof']).open('x') as h:json.dump(result,h,indent=2);h.write('\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
