#!/usr/bin/env python3
"""Rebuild endpoint memberships and verify family components with sparse graph traversal."""
import argparse,csv,gzip,json,hashlib,time,subprocess
from pathlib import Path
from collections import defaultdict,Counter
import psutil,numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--launch',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan);dep=json.loads(a.launch.read_text());assert ph==dep['plan_sha256']
    while True:
        try:
            p=psutil.Process(dep['pid'])
            if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines());assert state=={'ActiveState':'inactive','ExecMainStatus':'0'}
    for p,h in plan['pins'].items():assert sha(p)==h,p
    root=Path(plan['output']);r=json.loads((root/'receipt.json').read_text());assert r['status']=='complete_selected_control_family_dependency_inventory_pending_readback' and r['plan_sha256']==ph
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    selected={'target':set(),'background':set()};count=0
    with gzip.open(Path(plan['selection'])/'selections.tsv.gz','rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            selected['target'].add(row['target_id']);selected['background'].add(row['background_id']);count+=1
    expected=set();nodes=Counter();families=defaultdict(set)
    for role,ids in selected.items():
        seen=set()
        for line in (Path(plan['graph'])/(role+'_nodes.jsonl')).open():
            n=json.loads(line)
            if n['node_id'] not in ids:continue
            assert n['node_id'] not in seen;seen.add(n['node_id']);nodes[n['guide']+':'+role]+=1;families[n['guide']].add(n['family'])
            for end in ['a','b']:
                for kind,identity in [('gene',n['gene_'+end]),('model_version',json.dumps([n['model_id_'+end],n['version_'+end]],separators=(',',':'))),('sequence_sha256',n['sequence_sha256_'+end])]:
                    key=(n['guide'],role,n['node_id'],n['family'],end,kind,identity);assert key not in expected;expected.add(key)
        assert seen==ids
    total=len(expected);entities=defaultdict(set)
    columns=['guide','role','node_id','family','endpoint','entity_kind','entity_id']
    with gzip.open(root/'selected_endpoint_entities.tsv.gz','rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key=tuple(row[x] for x in columns);assert key in expected;expected.remove(key);entities[row['guide'],row['entity_kind'],row['entity_id']].add(row['family'])
    assert not expected and total==r['endpoint_entity_rows'] and dict(nodes)==r['selected_node_counts'] and count==r['selected_records'] and len(entities)==r['unique_entity_guide_keys']
    cross={}
    for row in csv.DictReader((root/'cross_family_entities.tsv').open(),delimiter='\t'):
        key=row['guide'],row['entity_kind'],row['entity_id'];assert key not in cross;cross[key]=set(json.loads(row['families']));assert len(cross[key])==int(row['family_count'])
    assert cross=={k:v for k,v in entities.items() if len(v)>1}
    actual={}
    for row in csv.DictReader((root/'family_components.tsv').open(),delimiter='\t'):
        key=row['guide'],row['entity_kind'],row['family'];assert key not in actual;actual[key]=(row['component_id'],int(row['component_families']))
    wanted={};summaries={}
    for guide in ['profile','mafft']:
        fs=sorted(families[guide]);indexes={x:i for i,x in enumerate(fs)}
        for kind in ['gene','model_version','sequence_sha256','combined']:
            relevant=[members for (g,k,e),members in entities.items() if g==guide and (kind=='combined' or kind==k)];ii=[];jj=[]
            for members in relevant:
                ix=sorted(indexes[x] for x in members)
                for j in ix[1:]:ii.append(ix[0]);jj.append(j)
            matrix=coo_matrix((np.ones(len(ii)),(ii,jj)),shape=(len(fs),len(fs))).tocsr();nc,labels=connected_components(matrix,directed=False)
            groups=defaultdict(list)
            for f,label in zip(fs,labels):groups[int(label)].append(f)
            for group in groups.values():
                for f in group:wanted[guide,kind,f]=(min(group),len(group))
            summaries[guide+':'+kind]=dict(families=len(fs),components=int(nc),multi_family_components=sum(len(g)>1 for g in groups.values()),maximum_component_families=max(map(len,groups.values()),default=0),cross_family_entities=sum(len(m)>1 for m in relevant))
    assert wanted==actual and summaries==r['component_summaries']
    for p,h in plan['pins'].items():assert sha(p)==h,p
    result=dict(status='passed_full_selected_control_family_dependency_readback',selected_records=count,endpoint_entity_rows=total,component_rows=len(actual),component_summaries=summaries,source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),scope='Every endpoint entity rederived from full selected-node sources; exact cross-family membership verified and every family component independently reconstructed by SciPy sparse connected-components traversal, not producer union-find. No phylogenetic independence or biological family equivalence claimed.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
