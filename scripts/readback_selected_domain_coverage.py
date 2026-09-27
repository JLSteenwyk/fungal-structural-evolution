#!/usr/bin/env python3
"""Independently rebuild all projected coverage summaries with NumPy and pandas."""
import argparse,csv,gzip,json,hashlib,itertools,time,subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import psutil


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--launch',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    plan=json.loads(a.plan.read_text());ph=sha(a.plan);dep=json.loads(a.launch.read_text());lh=sha(a.launch);assert dep['plan_sha256']==ph
    while True:
        try:
            p=psutil.Process(dep['pid'])
            if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines());assert state=={'ActiveState':'inactive','ExecMainStatus':'0'}
    assert sha(a.plan)==ph and sha(a.launch)==lh
    for p,h in plan['pins'].items():assert sha(p)==h,p
    out=Path(plan['output']);r=json.loads((out/'receipt.json').read_text());assert r['status']=='complete_selected_domain_coverage_projection_pending_independent_readback' and r['plan_sha256']==ph
    for name,h in r['artifacts'].items():assert sha(out/name)==h
    for name in ['inventory','qualification']:
        root=Path(plan[name]);rr=json.loads((root/'receipt.json').read_text());assert sha(root/'receipt.json')==r['source_bindings'][name]['receipt_sha256'] and sha(plan[name+'_audit'])==r['source_bindings'][name]['audit_sha256']
        for f,h in rr['artifacts'].items():assert sha(root/f)==h
    cfg={}
    for line in (Path(plan['inventory'])/'domain_configurations.jsonl').open():
        c=json.loads(line);assert c['domain_config_id'] not in cfg;cfg[c['domain_config_id']]=len(cfg)
    conditions=list(itertools.product(['alignment','envelope'],['full','plddt70','both'],plan['screens']));columns={x:i for i,x in enumerate(conditions)}
    statuses=['all_shared_domains_pass','some_shared_domains_pass','no_shared_domains_pass','no_shared_domain_comparison','identical_model_requires_separate_handling'];status_index={s:i for i,s in enumerate(statuses)}
    eligible=np.zeros((len(cfg),36),dtype=np.uint32);status=np.zeros((len(cfg),36),dtype=np.uint8);seen=np.zeros((len(cfg),36),dtype=bool)
    with gzip.open(Path(plan['qualification'])/'configuration_eligibility.tsv.gz','rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            i=cfg[row['domain_config_id']];j=columns[row['boundary'],row['mask_cohort'],row['screen']];assert not seen[i,j];seen[i,j]=True
            n=int(row['eligible_domain_count']);assert 0<=n<2**32;eligible[i,j]=n;status[i,j]=status_index[row['eligibility_status']]
    assert seen.all()
    fields=['target_id','background_id','policy','scenario_id','domain_config_id']
    selected=pd.read_csv(Path(plan['inventory'])/'selection_domain_links.tsv.gz',sep='\t',usecols=fields,dtype=str,keep_default_na=False)
    nodes=[]
    for line in Path(plan['target_nodes']).open():
        x=json.loads(line);nodes.append({k:x[k] for k in ['node_id','guide','family','taxon_id']})
    selected=selected.merge(pd.DataFrame(nodes),left_on='target_id',right_on='node_id',how='left',validate='many_to_one');assert len(selected)==2786912 and selected.notna().all().all()
    assert not selected.duplicated(['target_id','policy','scenario_id']).any()
    selected['config_index']=selected.domain_config_id.map(cfg);assert selected.config_index.notna().all()
    for col in ['background_id','family','taxon_id']:selected[col+'_code']=pd.factorize(selected[col],sort=True)[0]
    groups={key:group for key,group in selected.groupby(['guide','policy','scenario_id'],sort=False)}
    actual={}
    for row in csv.DictReader((out/'selected_domain_coverage.tsv').open(),delimiter='\t'):
        key=tuple(row[x] for x in ['guide','policy','scenario_id','boundary','mask_cohort','screen']);assert key not in actual;actual[key]=row
    scenarios=[x['scenario_id'] for x in json.loads(Path(plan['scenarios']).read_text())];checked=set();numeric=0
    for base in itertools.product(['profile','mafft'],plan['policies'],scenarios):
        group=groups.get(base);indices=np.array([],dtype=int) if group is None else group.config_index.to_numpy(dtype=int)
        code={x:np.array([],dtype=int) if group is None else group[x+'_code'].to_numpy(dtype=int) for x in ['taxon_id','family','background_id']}
        for condition,j in columns.items():
            key=base+condition;row=actual[key];checked.add(key);n=eligible[indices,j];ok=n>0;st=status[indices,j]
            expected=dict(selected_records=len(indices),usable_selected_records=int(ok.sum()),eligible_domain_occurrences=int(n.sum()))
            for source,name in [('taxon_id','taxa'),('family','families'),('background_id','backgrounds')]:
                expected['selected_'+name]=len(np.unique(code[source]));expected['usable_'+name]=len(np.unique(code[source][ok]))
            _,reuse=np.unique(code['background_id'][ok],return_counts=True);expected['maximum_usable_background_reuse']=int(reuse.max()) if len(reuse) else 0
            expected.update({s:int((st==i).sum()) for i,s in enumerate(statuses)})
            for col,value in expected.items():assert int(row[col])==value,(key,col,row[col],value);numeric+=1
    assert checked==set(actual) and len(checked)==r['summary_cells']==15552
    assert r['selected_records']==len(selected) and r['expanded_record_cells']==len(selected)*36
    for p,h in plan['pins'].items():assert sha(p)==h,p
    result=dict(status='passed_full_selected_domain_coverage_projection_readback',selected_records=len(selected),summary_cells=len(checked),numeric_values_checked=numeric,source_receipt_sha256=sha(out/'receipt.json'),script_sha256=sha(__file__),scope='Every summary independently rebuilt from full selected records, target metadata and audited qualification arrays, without SQL or producer helpers. Includes empty strata, every status, accession-occurrence count, taxon/family/background diversity and maximum control reuse. Not an independent-sample test or effect inference.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
