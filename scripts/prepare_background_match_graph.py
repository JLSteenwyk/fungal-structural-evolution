#!/usr/bin/env python3
"""Enumerate every admissible conserved-architecture control edge before choosing matches."""
import argparse,bisect,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from assess_architecture_matched_support import pair_signature,POLICIES
from run_ortholog_pair_guide_comparison import sha

def identity(values):return hashlib.sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest()
def eligible_range(values,d):
    if not math.isfinite(d) or d<0:raise ValueError('Invalid sequence distance')
    return bisect.bisect_left(values,d/2),bisect.bisect_right(values,d*2)
def features(left,right):
    return {field+'_'+side:model[field] for side,model in [('a',left),('b',right)] for field in ['model_id','version','length','mean_ca_plddt','fraction_ca_plddt_below50','sequence_sha256','sha256']}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify();source=json.loads(Path(plan['architecture_plan']).read_text());ar=json.loads(Path(plan['architecture_receipt']).read_text());proof=json.loads(Path(plan['architecture_readback']).read_text())
    assert proof['status']=='passed_full_architecture_matched_support_readback' and proof['producer_receipt_sha256']==sha(plan['architecture_receipt'])
    annotations={}
    for path in source['annotation_tables']:
        for line in open(path):
            r=json.loads(line);key=r['model_id'],r['version'];v={p:(tuple((h['pfam_accession'],h['pfam_type']) for h in hits),bool(hits) and all(h['conservative_architecture']==1 for h in hits)) for p,hits in r['policies'].items()}
            if key in annotations:assert annotations[key]==v
            annotations[key]=v
    models={}
    for path in plan['model_tables']:
        for line in open(path):
            m=json.loads(line);key=m['model_id'],m['version']
            if key in models:assert models[key]==m
            models[key]=m
    targets={}
    for r in csv.DictReader(open(source['target_links']),delimiter='\t'):
        key=r['guide'],r['gene_a'],r['gene_b'];assert key not in targets;targets[key]=r
    out=Path(plan['output']);out.mkdir(exist_ok=False);pools=defaultdict(list);backgrounds={}
    with (out/'background_nodes.jsonl').open('w') as f:
        for r in csv.DictReader(open(source['background_links']),delimiter='\t'):
            ends=[(r['model_id_'+s],int(r['version_'+s])) for s in ['a','b']]
            for guide in ['profile','mafft']:
                if r[guide+'_candidate_and_native_ortholog']!='1':continue
                ident=identity([guide,r['gene_a'],r['gene_b']]);assert ident not in backgrounds
                node=dict(node_id=ident,guide=guide,family=r[guide+'_family'],gene_a=r['gene_a'],gene_b=r['gene_b'],taxon_a=r['taxon_a'],taxon_b=r['taxon_b'],sequence_distance=float(r[guide+'_sequence_pair_distance']),pair_key=r['pair_key'],same_model=int(ends[0]==ends[1]),both_guides=int(r['candidate_and_native_ortholog_both']),both_unreported_parents=int(r['both_guides_and_parents_unreported']),**features(*[models[k] for k in ends]))
                assert math.isfinite(node['sequence_distance']) and node['sequence_distance']>=0
                backgrounds[ident]=node;f.write(json.dumps(node,separators=(',',':'))+'\n')
                for policy in POLICIES:
                    _,sig=pair_signature(*[annotations[k][policy] for k in ends])
                    if sig is not None:pools[guide,policy,node['family'],sig].append((node['sequence_distance'],ident))
    for v in pools.values():v.sort()
    distances={k:[x[0] for x in v] for k,v in pools.items()};seen=set();edge_counts=Counter();support_counts=Counter();nodes=set()
    with open(plan['architecture_table']) as f,(out/'edges.tsv').open('w') as e,(out/'target_nodes.jsonl').open('w') as t,(out/'target_policy_dispositions.tsv').open('w') as disposition:
        ew=csv.writer(e,delimiter='\t',lineterminator='\n');ew.writerow(['target_id','background_id','policy','focal_taxon','within_factor_1_25','within_factor_1_5'])
        dw=csv.writer(disposition,delimiter='\t',lineterminator='\n');dw.writerow(['target_id','policy','architecture_status','eligible_edges_within_factor_2'])
        for r in csv.DictReader(f,delimiter='\t'):
            if r['background_set']!='guide_native_ortholog':continue
            key=r['guide'],r['gene_a'],r['gene_b'];target=targets[key];ident=identity(key);policy=r['policy'];assert (ident,policy) not in seen;seen.add((ident,policy))
            d=float(r['target_sequence_distance']);ends=[(target['model_'+s],int(target['version_'+s])) for s in ['a','b']]
            if ident not in nodes:
                node=dict(node_id=ident,guide=key[0],family=target['family'],gene_a=key[1],gene_b=key[2],taxon_id=target['taxon_id'],gene_node=target['gene_node'],sequence_distance=d,pair_key=target['pair_key'],same_model=int(target['same_model']),**features(*[models[k] for k in ends]));t.write(json.dumps(node,separators=(',',':'))+'\n');nodes.add(ident)
            status,sig=pair_signature(*[annotations[k][policy] for k in ends]);assert status==r['target_architecture_status']
            poolkey=(key[0],policy,target['family'],sig);values=distances.get(poolkey,[]);lo,hi=eligible_range(values,d);assert hi-lo==int(r['within_factor_2_0']);counts=Counter()
            for distance,bid in pools.get(poolkey,[])[lo:hi]:
                b=backgrounds[bid];focal=int(target['taxon_id'] in [b['taxon_a'],b['taxon_b']]);flags=[int(d/factor<=distance<=d*factor) for factor in [1.25,1.5]]
                ew.writerow([ident,bid,policy,focal,*flags]);edge_counts[key[0]+'|'+policy]+=1
                for label,flag in [('1_25',flags[0]),('1_5',flags[1]),('2_0',1)]:counts[label]+=flag;counts['focal_'+label]+=flag*focal
            for label in ['1_25','1_5','2_0']:
                assert counts[label]==int(r['within_factor_'+label]) and counts['focal_'+label]==int(r['focal_within_factor_'+label])
            dw.writerow([ident,policy,status,hi-lo]);support_counts[status]+=1
    assert len(nodes)==len(targets)==ar['targets'] and len(seen)==4*len(targets) and sum(edge_counts.values())==plan['resources']['expected_edges']
    verify();r=dict(status='complete_background_match_graph_pending_readback',plan_sha256=ph,target_nodes=len(nodes),background_nodes=len(backgrounds),target_policy_dispositions=len(seen),edges=sum(edge_counts.values()),edge_counts=dict(edge_counts),target_architecture_counts=dict(support_counts),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All admissible edges within factor2 sequence distance for each guide and all four conservative shared-architecture policies. Exact factor1.25/1.5 and focal flags retained; background nodes retain both-guide and parent qualification flags. All targets and unsupported policy dispositions retained. Endpoint lengths and confidence preserved, not yet balanced or used to select matches. Reused backgrounds, genes, taxa and families are dependent; this is a candidate graph, not selected controls, independent replicates or effects. Structural outcomes not used in construction.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:r[k] for k in ['status','target_nodes','background_nodes','edges']}),flush=True)
if __name__=='__main__':main()
