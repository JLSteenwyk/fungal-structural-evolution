#!/usr/bin/env python3
"""Verify all graph nodes/edges and prove completeness against audited support counts."""
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha

def ident(values):return hashlib.sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest()
def readnodes(path):
    out={}
    for line in open(path):
        row=json.loads(line);key=row['node_id'];assert key not in out;out[key]=row
    return out

def architecture(node,policy,annotations):
    hits=[annotations[node['model_id_'+s],node['version_'+s]][policy] for s in ['a','b']]
    signatures=[tuple((h['pfam_accession'],h['pfam_type']) for h in part) for part in hits]
    if not hits[0] and not hits[1]:return 'neither_annotated',None
    if not hits[0] or not hits[1]:return 'one_unannotated',None
    if signatures[0]!=signatures[1]:return 'different_ordered_annotations',None
    if any(h['conservative_architecture']!=1 for part in hits for h in part):return 'shared_but_nonconservative',None
    return 'shared_conservative_architecture',signatures[0]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify();sp=json.loads(Path(plan['source_plan']).read_text());source=json.loads(Path(sp['architecture_plan']).read_text());root=Path(sp['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    assert r['status']=='complete_background_match_graph_pending_readback' and r['plan_sha256']==sha(plan['source_plan'])
    for path,h in sp['pins'].items():assert sha(path)==h,path
    for name,h in r['artifacts'].items():assert sha(root/name)==h,name
    proof=json.loads(Path(sp['architecture_readback']).read_text());assert proof['status']=='passed_full_architecture_matched_support_readback' and proof['producer_receipt_sha256']==sha(sp['architecture_receipt'])
    targets=readnodes(root/'target_nodes.jsonl');backgrounds=readnodes(root/'background_nodes.jsonl');models={}
    for path in sp['model_tables']:
        for line in open(path):
            m=json.loads(line);key=m['model_id'],m['version']
            if key in models:assert models[key]==m
            models[key]=m
    def addfeatures(node,ends):
        for side,key in zip(['a','b'],ends):
            m=models[key]
            for field in ['model_id','version','length','mean_ca_plddt','fraction_ca_plddt_below50','sequence_sha256','sha256']:node[field+'_'+side]=m[field]
        return node
    distances={}
    for row in csv.DictReader(open(source['support_table']),delimiter='\t'):
        if row['background_set']=='guide_native_ortholog':
            key=ident([row['guide'],row['gene_a'],row['gene_b']]);assert key not in distances;distances[key]=float(row['target_sequence_distance'])
    seen=set()
    for row in csv.DictReader(open(source['target_links']),delimiter='\t'):
        key=ident([row['guide'],row['gene_a'],row['gene_b']]);assert key not in seen;seen.add(key)
        ends=[(row['model_'+s],int(row['version_'+s])) for s in ['a','b']]
        expected=dict(node_id=key,guide=row['guide'],family=row['family'],gene_a=row['gene_a'],gene_b=row['gene_b'],taxon_id=row['taxon_id'],gene_node=row['gene_node'],sequence_distance=distances[key],pair_key=row['pair_key'],same_model=int(ends[0]==ends[1]))
        assert targets[key]==addfeatures(expected,ends)
    assert seen==set(targets)==set(distances)
    seen=set()
    for row in csv.DictReader(open(source['background_links']),delimiter='\t'):
        ends=[(row['model_id_'+s],int(row['version_'+s])) for s in ['a','b']]
        for guide in ['profile','mafft']:
            if row[guide+'_candidate_and_native_ortholog']!='1':continue
            key=ident([guide,row['gene_a'],row['gene_b']]);assert key not in seen;seen.add(key)
            expected=dict(node_id=key,guide=guide,family=row[guide+'_family'],gene_a=row['gene_a'],gene_b=row['gene_b'],taxon_a=row['taxon_a'],taxon_b=row['taxon_b'],sequence_distance=float(row[guide+'_sequence_pair_distance']),pair_key=row['pair_key'],same_model=int(ends[0]==ends[1]),both_guides=int(row['candidate_and_native_ortholog_both']),both_unreported_parents=int(row['both_guides_and_parents_unreported']))
            assert backgrounds[key]==addfeatures(expected,ends)
    assert seen==set(backgrounds);del models,distances
    annotations={}
    for path in source['annotation_tables']:
        for line in open(path):
            row=json.loads(line);key=row['model_id'],row['version'];v=row['policies']
            if key in annotations:assert annotations[key]==v
            annotations[key]=v
    policies=['alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore']
    signatures={};statuses={}
    for nodes in [targets,backgrounds]:
        for key,node in nodes.items():
            for policy in policies:statuses[key,policy],signatures[key,policy]=architecture(node,policy,annotations)
    del annotations
    edges=set();counts=defaultdict(Counter);ec=Counter()
    for row in csv.DictReader(open(root/'edges.tsv'),delimiter='\t'):
        key=row['target_id'],row['background_id'],row['policy'];assert key not in edges;edges.add(key);tid,bid,policy=key;t=targets[tid];b=backgrounds[bid]
        assert policy in policies and t['guide']==b['guide'] and t['family']==b['family']
        assert signatures[tid,policy] is not None and signatures[tid,policy]==signatures[bid,policy]
        d=t['sequence_distance'];db=b['sequence_distance'];assert d/2<=db<=2*d
        focal=int(t['taxon_id'] in [b['taxon_a'],b['taxon_b']]);assert int(row['focal_taxon'])==focal
        flags={label:int(d/factor<=db<=d*factor) for label,factor in [('1_25',1.25),('1_5',1.5),('2_0',2)]}
        assert all(int(row['within_factor_'+label])==flags[label] for label in ['1_25','1_5'])
        sets=['guide_native_ortholog']
        if b['both_guides']:sets.append('both_guides_native_ortholog')
        if b['both_unreported_parents']:sets.append('both_guides_unreported_parents')
        for setting in sets:
            c=counts[tid,policy,setting]
            for label,flag in flags.items():c[label]+=flag;c['focal_'+label]+=flag*focal
        ec[t['guide']+'|'+policy]+=1
    assert len(edges)==r['edges'] and dict(ec)==r['edge_counts'];del edges
    checked=set();statuscounts=Counter()
    for row in csv.DictReader(open(root/'target_policy_dispositions.tsv'),delimiter='\t'):
        key=row['target_id'],row['policy'];assert key not in checked;checked.add(key)
        assert row['architecture_status']==statuses[key] and int(row['eligible_edges_within_factor_2'])==counts.get((*key,'guide_native_ortholog'),{}).get('2_0',0);statuscounts[row['architecture_status']]+=1
    assert checked=={(t,p) for t in targets for p in policies} and dict(statuscounts)==r['target_architecture_counts']
    checked=set()
    for row in csv.DictReader(open(sp['architecture_table']),delimiter='\t'):
        tid=ident([row['guide'],row['gene_a'],row['gene_b']]);key=tid,row['policy'],row['background_set'];assert key not in checked;checked.add(key);c=counts.get(key,{})
        assert statuses[key[:2]]==row['target_architecture_status']
        for label in ['1_25','1_5','2_0']:
            assert c.get(label,0)==int(row['within_factor_'+label]) and c.get('focal_'+label,0)==int(row['focal_within_factor_'+label])
    assert len(checked)==12*len(targets) and len(targets)==r['target_nodes'] and len(backgrounds)==r['background_nodes'] and 4*len(targets)==r['target_policy_dispositions']
    verify();assert sha(rp)==rh
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    result=dict(status='passed_full_background_match_graph_readback',plan_sha256=ph,producer_receipt_sha256=rh,target_nodes=len(targets),background_nodes=len(backgrounds),edges=r['edges'],target_policy_dispositions=r['target_policy_dispositions'],support_rows_checked=len(checked),scope='Every target/background field rejoined to source genes, model descriptors and guide distances. All unique edges directly checked for guide/family/conservative architecture/range eligibility and flags; exact counts reconciled to every audited support row in all three qualification sets and three distance ranges, proving completeness. All unsupported target/policy records verified. No matching selection, balance or effects established.')
    Path(plan['output']).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
