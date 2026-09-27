#!/usr/bin/env python3
"""Audit cross-family gene, model and sequence reuse among every selected control."""
import argparse,csv,gzip,json,itertools
from pathlib import Path
from collections import defaultdict,Counter
from screen_duplication_domain_alignment_coverage import sha


def components(families,memberships):
    parent={f:f for f in families}
    def find(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    for members in memberships:
        members=sorted(members)
        for f in members[1:]:
            a,b=find(members[0]),find(f)
            if a!=b:parent[max(a,b)]=min(a,b)
    groups=defaultdict(list)
    for f in sorted(families):groups[find(f)].append(f)
    return [v for k,v in sorted(groups.items())]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();selection=Path(plan['selection']);sr=json.loads((selection/'receipt.json').read_text());proof=json.loads(Path(plan['selection_audit']).read_text());assert proof['producer_receipt_sha256']==sha(selection/'receipt.json') and proof['status']=='passed_full_background_control_selection_readback'
    assert sha(selection/'selections.tsv.gz')==sr['artifacts']['selections.tsv.gz']
    selected={'target':set(),'background':set()};rows=0
    with gzip.open(selection/'selections.tsv.gz','rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            for role in selected:selected[role].add(row[role+'_id'])
            rows+=1
    assert rows==sr['selected_records'];graph=Path(plan['graph']);gr=json.loads((graph/'receipt.json').read_text());ga=json.loads(Path(plan['graph_audit']).read_text());assert ga['producer_receipt_sha256']==sha(graph/'receipt.json')
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);entities=defaultdict(set);families=defaultdict(set);node_counts=Counter();endpoint_count=0
    with gzip.open(out/'selected_endpoint_entities.tsv.gz','wt') as f:
        writer=csv.DictWriter(f,['guide','role','node_id','family','endpoint','entity_kind','entity_id'],delimiter='\t',lineterminator='\n');writer.writeheader()
        for role,ids in selected.items():
            path=graph/(role+'_nodes.jsonl');assert sha(path)==gr['artifacts'][path.name];seen=set()
            for line in path.open():
                n=json.loads(line)
                if n['node_id'] not in ids:continue
                assert n['node_id'] not in seen;seen.add(n['node_id']);guide=n['guide'];family=n['family'];families[guide].add(family);node_counts[guide+':'+role]+=1
                for end in ['a','b']:
                    identity={'gene':n['gene_'+end],'model_version':json.dumps([n['model_id_'+end],n['version_'+end]],separators=(',',':')),'sequence_sha256':n['sequence_sha256_'+end]}
                    for kind,value in identity.items():
                        entities[guide,kind,value].add(family);writer.writerow(dict(guide=guide,role=role,node_id=n['node_id'],family=family,endpoint=end,entity_kind=kind,entity_id=value));endpoint_count+=1
            assert seen==ids
    with (out/'cross_family_entities.tsv').open('w') as f:
        writer=csv.DictWriter(f,['guide','entity_kind','entity_id','families','family_count'],delimiter='\t',lineterminator='\n');writer.writeheader()
        for (guide,kind,identity),fs in sorted(entities.items()):
            if len(fs)>1:writer.writerow(dict(guide=guide,entity_kind=kind,entity_id=identity,families=json.dumps(sorted(fs),separators=(',',':')),family_count=len(fs)))
    summaries={};output=[]
    for guide in ['profile','mafft']:
        for kind in ['gene','model_version','sequence_sha256','combined']:
            maps=[fs for (g,k,e),fs in entities.items() if g==guide and (kind=='combined' or k==kind)]
            groups=components(families[guide],maps)
            for fs in groups:
                for family in fs:output.append(dict(guide=guide,entity_kind=kind,family=family,component_id=fs[0],component_families=len(fs)))
            summaries[guide+':'+kind]=dict(families=len(families[guide]),components=len(groups),multi_family_components=sum(len(x)>1 for x in groups),maximum_component_families=max(map(len,groups),default=0),cross_family_entities=sum(len(x)>1 for x in maps))
    with (out/'family_components.tsv').open('w') as f:
        writer=csv.DictWriter(f,list(output[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(output)
    verify();result=dict(status='complete_selected_control_family_dependency_inventory_pending_readback',plan_sha256=ph,selected_records=rows,selected_node_counts=dict(node_counts),endpoint_entity_rows=endpoint_count,unique_entity_guide_keys=len(entities),component_summaries=summaries,artifacts={p.name:sha(p) for p in out.iterdir()},scope='All unique nodes participating in any selected scenario/policy, before structural qualification. Shared genes, exact model/version identities and sequence hashes define separate and combined family-connection graphs within each guide. Repeated policies/scenarios do not multiply node membership. No claim that connected families are biologically the same, that components are independent phylogenetically, or that this is an optimal resampling scheme.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
