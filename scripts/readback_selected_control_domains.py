#!/usr/bin/env python3
"""Verify every selected-control link and independently reconstruct domain configurations."""
import argparse,csv,gzip,json,hashlib,itertools,time,subprocess
from pathlib import Path
from collections import Counter
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
    for p,h in plan['pins'].items():assert sha(p)==h
    root=Path(plan['output']);r=json.loads((root/'receipt.json').read_text());assert r['status']=='complete_selected_control_domain_inventory_pending_readback' and r['plan_sha256']==ph
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    indexes={}
    for label,spec in plan['domains'].items():
        index={}
        for x in csv.DictReader((Path(spec['root'])/'domain_pair_links.tsv').open(),delimiter='\t'):
            if label=='target' and x['scope'] not in ['primary','primary+reference']:continue
            key=x['pair_key'],x['policy'];sub=index.setdefault(key,{})
            slot=x['pfam_accession'],x['boundary'];assert slot not in sub;sub[slot]=x['domain_pair_key']
        indexes[label]=index
    configs={}
    for line in (root/'domain_configurations.jsonl').open():
        c=json.loads(line);cid=c['domain_config_id'];assert cid not in configs;configs[cid]=c
        identity=[c['target_pair_key'],c['background_pair_key'],c['policy']]
        assert cid==hashlib.sha256(json.dumps(identity,separators=(',',':')).encode()).hexdigest()
        maps=[indexes[label].get((c[label+'_pair_key'],c['policy']),{}) for label in ['target','background']]
        families=[{x for x,y in m} for m in maps]
        assert c['target_accessions']==sorted(families[0]) and c['background_accessions']==sorted(families[1])
        shared=families[0].intersection(families[1]);assert c['shared_accessions']==sorted(shared)
        expected={(pf,bound,maps[0][pf,bound],maps[1][pf,bound]) for pf in shared for bound in ['alignment','envelope']}
        actual={(x['pfam_accession'],x['boundary'],x['target_domain_pair'],x['background_domain_pair']) for x in c['matches']}
        assert actual==expected and len(actual)==len(c['matches'])
        status='identical_model_requires_separate_handling' if c['target_same_model'] or c['background_same_model'] else 'shared_domain_comparisons' if expected else 'no_shared_domain_comparison'
        assert c['status']==status
    nodes={}
    for label in ['target','background']:
        nodes[label]={}
        for line in (Path(plan['graph'])/(label+'_nodes.jsonl')).open():
            n=json.loads(line);assert n['node_id'] not in nodes[label];nodes[label][n['node_id']]={k:n[k] for k in ['pair_key','same_model','guide','family']}
    count=0;used=set();totals=Counter();seen=set()
    with gzip.open(Path(plan['selection'])/'selections.tsv.gz','rt') as f,gzip.open(root/'selection_domain_links.tsv.gz','rt') as g:
        for source,actual in itertools.zip_longest(csv.DictReader(f,delimiter='\t'),csv.DictReader(g,delimiter='\t')):
            assert source is not None and actual is not None
            cid=actual.pop('domain_config_id');assert actual==source
            key=source['target_id'],source['policy'],source['scenario_id'];assert key not in seen;seen.add(key)
            c=configs[cid];used.add(cid)
            for label in ['target','background']:
                n=nodes[label][source[label+'_id']]
                assert c[label+'_pair_key']==n['pair_key'] and c[label+'_same_model']==n['same_model']
            t=nodes['target'][source['target_id']];b=nodes['background'][source['background_id']];assert t['guide']==b['guide'] and t['family']==b['family'] and c['policy']==source['policy']
            totals[t['guide']+'|'+source['policy']+'|'+source['scenario_id']+'|'+c['status']]+=1;count+=1
    assert used==set(configs) and len(configs)==r['distinct_configurations']
    assert count==r['selected_records'] and dict(totals)==r['selected_status_counts']
    assert dict(Counter(c['status'] for c in configs.values()))==r['configuration_status_counts']
    for p,h in plan['pins'].items():assert sha(p)==h
    result=dict(status='passed_full_selected_control_domain_inventory_readback',selected_records=count,distinct_configurations=len(configs),configuration_status_counts=r['configuration_status_counts'],source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),scope='Every selected source field and node/policy binding checked; complete shared Pfam/boundary pair sets reconstructed without producer helpers; every configuration used and all aggregate counts verified. Not coverage, structural effect or independent biological validation.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
