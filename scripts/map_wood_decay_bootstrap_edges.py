#!/usr/bin/env python3
"""Map all optimal ecological edge assignments across the frozen bootstrap ensemble."""
import argparse,csv,gzip,json,hashlib,fcntl
from collections import Counter,defaultdict
from pathlib import Path
from Bio import Phylo
from ecology_optimal_edge_states import edge_costs

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text());ph=sha(a.plan)
 def verify():
  assert sha(a.plan)==ph
  for path,h in p['pins'].items():assert sha(path)==h,path
 verify();out=Path(p['output']);out.mkdir(parents=True,exist_ok=True);lock=(out/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 with open(p['evidence']) as f:evidence=list(csv.DictReader(f,delimiter='\t'))
 states={r['taxon_id']:int(r['state']=='brown_rot') for r in evidence if r['state'] in ['white_rot','brown_rot']};assert len(states)==6
 with open(p['manifest']) as f:taxa=sorted(r['taxon_id'] for r in csv.DictReader(f,delimiter='\t'))
 assert len(taxa)==len(set(taxa));alltaxa=set(taxa);splits={};aggregate=defaultdict(Counter);summaries=[];shards={};total=0
 for source in p['trees']:
  sourcehash=sha(source['bootstrap']);r=json.loads(Path(source['receipt']).read_text());proof=json.loads(Path(source['readback']).read_text())
  assert proof['source_receipt_sha256']==sha(source['receipt']) and proof['bootstrap_trees']==1000
  assert sourcehash==r['artifacts'][Path(source['bootstrap']).name]
  for number,tree in enumerate(Phylo.parse(source['bootstrap'],'newick'),1):
   tips=[n.name for n in tree.get_terminals()];assert len(tips)==len(taxa) and set(tips)==alltaxa
   nodes=list(tree.find_clades(order='preorder'));index={n:i for i,n in enumerate(nodes)};desc={}
   for node in reversed(nodes):desc[node]={node.name} if node.is_terminal() else set().union(*(desc[c] for c in node.clades))
   identifiers={}
   for node in nodes[1:]:
    side=sorted(desc[node]);other=sorted(alltaxa-desc[node]);canonical=min([side,other],key=lambda x:(len(x),x));key=hashlib.sha256(json.dumps(canonical,separators=(',',':')).encode()).hexdigest();identifiers[node]=key
    if key in splits:assert splits[key]==canonical
    else:splits[key]=canonical
   assert len(set(identifiers.values()))==len(nodes)-1
   shard=out/f"{source['label']}-{number:04d}.json.gz";sr=Path(str(shard)+'.receipt.json')
   if sr.exists():
    saved=json.loads(sr.read_text());assert saved['plan_sha256']==ph and saved['source_sha256']==sourcehash and saved['index']==number and saved['artifact_sha256']==sha(shard)
    with gzip.open(shard,'rt') as f:data=json.load(f)
   else:
    assert not shard.exists();data=[]
    for scenario in p['scenarios']:
     coding=scenario['label'];current=dict(states,**scenario['uncertain_assignments']);score,edges=edge_costs(tree,current);inf=len(nodes)+1
     entries=[]
     for parent,child,costs in edges:
      same=any(costs[i]==score for i in [0,3]);change=any(costs[i]==score for i in [1,2]);status='optional_change' if same and change else 'required_change' if change else 'no_change_in_any_optimum'
      entries.append([identifiers[child],index[parent],index[child],status,[x if x<inf else None for x in costs]])
     data.append(dict(coding=coding,minimum_changes=score,edges=entries))
    partial=Path(str(shard)+'.partial')
    with gzip.open(partial,'wt') as f:json.dump(data,f,separators=(',',':'))
    partial.replace(shard);sr.write_text(json.dumps(dict(plan_sha256=ph,source_sha256=sourcehash,index=number,artifact_sha256=sha(shard)))+'\n')
   assert {d['coding'] for d in data}=={x['label'] for x in p['scenarios']}
   for d in data:
    counts=Counter();assert len(d['edges'])==len(nodes)-1
    for key,parent,child,status,costs in d['edges']:
     assert identifiers[nodes[child]]==key and nodes[child] in nodes[parent].clades
     counts[status]+=1;aggregate[source['label'],d['coding'],key]['present']+=1;aggregate[source['label'],d['coding'],key][status]+=1;total+=1
    summaries.append(dict(tree_source=source['label'],tree_index=number,coding=d['coding'],minimum_changes=d['minimum_changes'],edges=len(d['edges']),**{k:counts[k] for k in ['required_change','optional_change','no_change_in_any_optimum']}))
   shards[shard.name]=sha(shard)
   if number%25==0:
    state=dict(source=source['label'],completed_trees=number,total_source_trees=1000,edge_rows=total);(out/'state.json').write_text(json.dumps(state)+'\n');print(json.dumps(state),flush=True)
  assert number==1000
 def write(name,rows):
  with (out/name).open('w') as f:
   w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 write('tree_summaries.tsv',summaries)
 output=[]
 for (source,coding,split),counts in sorted(aggregate.items()):
  assert sum(counts[k] for k in ['required_change','optional_change','no_change_in_any_optimum'])==counts['present']<=1000
  output.append(dict(tree_source=source,coding=coding,split_id=split,canonical_side=';'.join(splits[split]),ensemble_trees=1000,**{k:counts[k] for k in ['present','required_change','optional_change','no_change_in_any_optimum']}))
 write('split_frequencies.tsv',output);verify()
 result=dict(status='complete_wood_decay_bootstrap_edge_mapping_pending_readback',plan_sha256=ph,trees=2000,coding_tree_combinations=len(summaries),edge_rows=total,split_summary_rows=len(output),shards=shards,artifacts={name:sha(out/name) for name in ['tree_summaries.tsv','split_frequencies.tsv']},scope='All 2000 saved ultrafast bootstrap trees, five explicit decay coding scenarios. Per-edge constrained costs retained in hashed shards. Split presence and required/optional counts have separate denominators; frequencies are conditional topology sensitivity, not posterior transition probabilities, rooted origins, independent replication or trait/model uncertainty. Independent full bootstrap readback remains pending.')
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print('Complete',len(summaries),total,flush=True)
if __name__=='__main__':main()
