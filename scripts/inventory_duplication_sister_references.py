#!/usr/bin/env python3
"""Inventory immediate sister-clade model references without asserting biological outgroups."""
import argparse,csv,json,math,sqlite3
from collections import Counter,defaultdict
from io import StringIO
from pathlib import Path
from Bio import Phylo
from run_ortholog_pair_guide_comparison import sha


def tree_index(newick):
    tree=Phylo.read(StringIO(newick),'newick');nodes={};parents={};distance={};stack=[(tree.root,None,0.)]
    while stack:
        node,parent,d=stack.pop()
        if not node.name or node.name in nodes:raise ValueError('Ambiguous node')
        if parent is not None:
            if node.branch_length is None or not math.isfinite(node.branch_length) or node.branch_length<0:raise ValueError('Invalid branch')
            d+=node.branch_length
        nodes[node.name]=node;parents[node.name]=parent;distance[node.name]=d
        stack.extend((child,node.name,d) for child in node.clades)
    return nodes,parents,distance


def reference_summary(row,index,duplicated_nodes,models,taxa):
    nodes,parents,distance=index;node=nodes[row['gene_node']];focal=row['taxon_id'];genes=sorted([c.name for c in node.clades if not c.clades])
    if len(node.clades)!=2 or genes!=sorted([row['gene_a'],row['gene_b']]):raise ValueError('Candidate no longer exact terminal pair')
    parent=parents[node.name]
    result=dict(parent_node=parent or '',parent_reported_duplication=int(parent in duplicated_nodes),parent_children=0,sister_genes=0,sister_taxa=0,sister_focal_taxon_genes=0,modeled_nonfocal_sister_genes=0,modeled_nonfocal_sister_taxa=0,nearest_reference_genes='[]',chosen_reference_gene='',reference_model='',reference_version='',reference_distance_from_duplicate_node='',distance_a_to_reference='',distance_b_to_reference='',duplicate_pair_sequence_distance=nodes[row['gene_a']].branch_length+nodes[row['gene_b']].branch_length,status='no_parent')
    if parent is None:return result
    pn=nodes[parent];result['parent_children']=len(pn.clades);stack=[(c,c.branch_length) for c in pn.clades if c is not node];staxa=set();mtaxa=set();available=[]
    while stack:
        tip,path_length=stack.pop()
        if tip.clades:stack.extend((c,path_length+c.branch_length) for c in tip.clades);continue
        taxon=tip.name.split('_',1)[0]
        if taxon not in taxa:raise ValueError('Unknown taxon prefix in resolved tree')
        result['sister_genes']+=1;staxa.add(taxon)
        if taxon==focal:result['sister_focal_taxon_genes']+=1
        elif tip.name in models:
            # Node and reference have this immediate parent as their MRCA.
            d=node.branch_length+path_length
            available.append((d,tip.name));mtaxa.add(taxon)
    result.update(sister_taxa=len(staxa),modeled_nonfocal_sister_genes=len(available),modeled_nonfocal_sister_taxa=len(mtaxa))
    if available:
        best=min(d for d,g in available);ties=sorted(g for d,g in available if math.isclose(d,best,rel_tol=0,abs_tol=1e-12));chosen=ties[0];m=models[chosen]
        actual=next(d for d,g in available if g==chosen)
        result.update(nearest_reference_genes=json.dumps(ties),chosen_reference_gene=chosen,reference_model=m[0],reference_version=m[1],reference_distance_from_duplicate_node=actual,distance_a_to_reference=actual+nodes[row['gene_a']].branch_length,distance_b_to_reference=actual+nodes[row['gene_b']].branch_length)
    if parent in duplicated_nodes:result['status']='parent_reported_duplication'
    elif len(pn.clades)!=2:result['status']='parent_not_bifurcating'
    elif result['sister_focal_taxon_genes']:result['status']='focal_taxon_in_sister_clade'
    elif not available:result['status']='no_modeled_nonfocal_sister'
    else:result['status']='provisional_reference_available'
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed pin '+path)
    verify();review=Path(plan['review']);rr=json.loads((review/'receipt.json').read_text());rp=json.loads(Path(plan['review_plan']).read_text())
    if rr['status']!='complete_duplication_candidate_join_and_tree_review' or rr['plan_sha256']!=sha(plan['review_plan']):raise ValueError('Candidate review mismatch')
    if str(plan['bridge'])!=rp['bridge'] or sha(plan['bridge'])!=rp['pins'][rp['bridge']]:raise ValueError('Bridge differs from reviewed source')
    db=sqlite3.connect('file:'+str(Path(plan['bridge']).resolve())+'?mode=ro',uri=True);models={t+'_'+p:(m,v) for t,p,m,v in db.execute('SELECT taxon_id,protein_id,model_id,version FROM structures')};db.close()
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);summaries=[]
    for entry in rp['guides']:
        guide=entry['guide'];path=review/(guide+'_candidate_tree_checks.tsv')
        if sha(path)!=rr['artifacts'][path.name]:raise ValueError('Changed reviewed candidates')
        with path.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        groups=defaultdict(list)
        for r in rows:
            if r['tree_status']!='exact_reported_pair':raise ValueError('Unresolved candidate tree')
            groups[r['family']].append(r)
        dup=defaultdict(set)
        with Path(entry['events']).open() as f:
            for r in csv.DictReader(f,delimiter='\t'):
                if r['Orthogroup'] in groups:dup[r['Orthogroup']].add(r['Gene Tree Node'])
        tree_audit=json.loads(Path(entry['tree_readback']).read_text());species=next(Path(p) for p in tree_audit['input_hashes'] if Path(p).name=='SpeciesIDs.txt')
        taxa={line.split(': ',1)[1].rsplit('.',1)[0] for line in species.read_text().splitlines()}
        if any('_' in t for t in taxa):raise ValueError('Taxon prefix needs explicit mapping')
        counts=Counter();seen=set();writer=None
        with Path(entry['trees']).open() as f,(out/(guide+'_sister_references.tsv')).open('w') as output:
            for line in f:
                family,newick=line.rstrip().split(': ',1)
                if family not in groups:continue
                if family in seen:raise ValueError('Repeated tree')
                seen.add(family);index=tree_index(newick)
                for r in groups[family]:
                    result=reference_summary(r,index,dup[family],models,taxa);counts[result['status']]+=1;data={**r,**result}
                    if writer is None:writer=csv.DictWriter(output,fieldnames=list(data),delimiter='\t',lineterminator='\n');writer.writeheader()
                    writer.writerow(data)
                if len(seen)%1000==0:print(guide,'families',len(seen),flush=True)
        if seen!=set(groups) or sum(counts.values())!=len(rows):raise ValueError('Incomplete candidate scope')
        summaries.append(dict(guide=guide,candidates=len(rows),families=len(seen),counts=dict(counts)));print(json.dumps(summaries[-1]),flush=True)
    verify();result=dict(status='complete_duplication_sister_reference_inventory',plan_sha256=ph,guides=summaries,artifacts={p.name:sha(p) for p in out.iterdir()},scope='Immediate sister-clade availability for every reviewed terminal duplicate pair. Provisional references require a bifurcating parent not reported as duplication and no focal-taxon genes in the sister clade. Nearest modeled nonfocal tip selected by fixed sequence-tree path length, ties within 1e-12 retained with lexical representative. These are extant reference candidates, not validated orthologous outgroups or ancestral states. Parent not reported as duplication does not establish speciation. Distances are model-relative sequence divergence, not time or structural asymmetry. Reference coordinate validation, domain/functional checks and guide/reference sensitivity remain pending.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    csv.field_size_limit(32*1024*1024);main()
