#!/usr/bin/env python3
"""Rebuild every fitting-input fingerprint through independent SQL joins."""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path
import duckdb
import numpy as np
import pandas as pd
import psutil
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser()
    for name in ['plan','launch','output']:ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args(); plan=json.loads(args.plan.read_text()); ph=sha(args.plan)
    launch=json.loads(args.launch.read_text());lh=sha(args.launch)
    assert launch['plan_sha256']==ph
    while True:
        try:
            process=psutil.Process(launch['pid'])
            if process.create_time()!=launch['created'] or process.status()==psutil.STATUS_ZOMBIE:break
            assert process.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines())
    assert state=={'ActiveState':'inactive','ExecMainStatus':'0'},state
    assert sha(args.plan)==ph and sha(args.launch)==lh
    for path,digest in plan['pins'].items():assert sha(path)==digest,path
    root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['plan_sha256']==ph
    for name,digest in receipt['artifacts'].items():assert sha(root/name)==digest,name
    mapping=pd.read_csv(root/'full_setting_fit_map.tsv',sep='\t',dtype={'target_order':str,'background_order':str})
    keys=['guide','policy','scenario_id'];settings=['boundary','mask','cohort','screen','target_order','background_order']
    assert len(mapping)==82944 and not mapping.duplicated(keys+settings).any()
    recipes={}
    for line in (root/'unique_fit_recipes.jsonl').open():
        row=json.loads(line);assert row['fit_input_id'] not in recipes
        recipes[row['fit_input_id']]=row
    source=Path(plan['summaries'])
    summary=pd.read_csv(source/'record_summary.tsv',sep='\t',dtype={'target_order':str,'background_order':str})
    pd.testing.assert_series_equal(mapping.set_index(keys+settings).records.sort_index(),
                                  summary.set_index(keys+settings).matched_records.rename('records').sort_index())
    db=duckdb.connect();db.execute('SET threads=1');db.execute("SET memory_limit='8GB'")
    db.execute("CREATE TABLE target AS SELECT node_id,guide,family_component FROM read_csv(?,delim='\t',header=true,all_varchar=true) WHERE role='target'",[plan['nodes']])
    db.execute("CREATE TABLE pairs AS SELECT p.target_id,p.background_id,t.guide,sha256(to_json([p.target_id,p.background_id,t.family_component,p.species_pattern_id])) row_identity FROM read_csv(?,delim='\t',header=true,all_varchar=true) p JOIN target t ON p.target_id=t.node_id",[plan['pairs']])
    db.execute("CREATE TABLE selections AS SELECT s.*,p.guide,p.row_identity FROM read_csv(?,delim='\t',header=true,all_varchar=true) s JOIN pairs p USING(target_id,background_id)",[plan['selections']])
    assert db.execute('SELECT count(*) FROM selections').fetchone()[0]==2786912
    numeric=['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference']
    checked=0;seen=set();recipe_seen=set()
    for part in json.loads((source/'partition_manifest.json').read_text()):
        path=source/part['path'];assert sha(path)==part['sha256']
        data=db.execute('SELECT s.guide,s.policy,s.scenario_id,s.target_id,s.row_identity,'+','.join('v.'+k for k in numeric)+' FROM selections s JOIN read_parquet(?) v USING(domain_config_id) ORDER BY s.guide,s.policy,s.scenario_id,s.target_id',[str(path)]).df()
        expected=mapping[mapping.partition_index.eq(part['index'])].set_index(keys)
        assert len(expected)==432
        local=0
        for key,frame in data.groupby(keys,sort=False):
            row=expected.loc[key]
            for name in settings:assert row[name]==str(part[name])
            assert len(frame)==int(row.records)
            matrix=np.ascontiguousarray(frame[numeric].to_numpy(),dtype='<f8');assert np.isfinite(matrix).all()
            matrix[matrix==0]=0.
            spec=dict(records=len(frame),numeric_columns=numeric,
                      values_sha256=hashlib.sha256(matrix.tobytes()).hexdigest(),
                      ordered_identity_sha256=hashlib.sha256(np.asarray(frame.row_identity,dtype='S64').tobytes()).hexdigest(),
                      model_specification=plan['model_specification'])
            fingerprint=hashlib.sha256(json.dumps(spec,sort_keys=True,separators=(',',':')).encode()).hexdigest()
            assert fingerprint==row.fit_input_id
            recipe=recipes[fingerprint]
            for name,value in spec.items():assert recipe[name]==value
            if recipe['partition_index']==part['index'] and all(recipe[k]==v for k,v in zip(keys,key)):
                recipe_seen.add(fingerprint)
            seen.add(fingerprint);local+=1;checked+=1
        assert local==432
        print('Independently checked fit inventory',part['index']+1,'/ 192',flush=True)
    assert checked==receipt['full_settings']==82944 and seen==recipe_seen==set(recipes)
    assert len(recipes)==receipt['unique_record_inputs']
    assert receipt['full_tree_fits']==checked*5 and receipt['unique_tree_fits']==len(recipes)*5
    assert receipt['maximum_unique_optimizer_attempts']==len(recipes)*5*22
    assert receipt['record_occurrences']==int(mapping.records.sum())
    assert receipt['unique_input_record_occurrences']==sum(r['records'] for r in recipes.values())
    for path,digest in plan['pins'].items():assert sha(path)==digest,path
    result=dict(status='passed_full_matched_fit_inventory_readback',settings=checked,unique_record_inputs=len(recipes),
                source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),
                scope='Independent SQL joins and ordered identity construction reproduce all input hashes, recipes, labels and denominators. No approximate equivalence, effect fit or inference.')
    with args.output.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
