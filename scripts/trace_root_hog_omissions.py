#!/usr/bin/env python3
"""Trace every unflagged root omission against emitted HOG parent clades."""
import argparse
import csv
from io import StringIO
import json
from pathlib import Path
from Bio import Phylo
from assess_small_family_output_exposure import sha


def trace_family(newick, hog_parents, flagged, targets):
    tree=Phylo.read(StringIO(newick),'newick');nodes={};parents={};counts={};hog_counts={};paths={};contents={k:set() for k in hog_parents}
    stack=[(tree.root,None)]
    while stack:
        node,active=stack.pop()
        if node.name in nodes:raise ValueError('Repeated tree node label')
        nodes[node.name]=node
        if node.name in hog_parents:
            if active is not None:raise ValueError('Nested root HOG parent clades')
            active=node.name
        if node.clades:
            for child in node.clades:parents[child.name]=node.name;stack.append((child,active))
        else:
            paths[node.name]=active
            if active is not None:contents[active].add(node.name)
    for parent,expected in hog_parents.items():
        if parent not in nodes or not nodes[parent].clades:raise ValueError('Missing HOG parent node')
        if contents[parent]-flagged!=expected:raise ValueError('HOG parent descendants differ from emitted membership')
    if not targets<=paths.keys():raise ValueError('Missing target gene in resolved tree')
    if targets.intersection(set().union(*hog_parents.values())):raise ValueError('Target already assigned to root HOG')
    # Explicit postorder traversal avoids recursion on giant families.
    stack=[(tree.root,False)]
    while stack:
        node,visited=stack.pop()
        if not visited:
            stack.append((node,True));stack.extend((ch,False) for ch in node.clades);continue
        counts[node.name]=sum(counts[ch.name] for ch in node.clades) if node.clades else 1
        hog_counts[node.name]=(node.name in hog_parents)+sum(hog_counts[ch.name] for ch in node.clades)
    rows=[]
    for target in sorted(targets):
        parent=parents.get(target);join=parent;distance=1
        while join is not None and not hog_counts[join]:join=parents.get(join);distance+=1
        rows.append(dict(gene_label=target,inside_emitted_parent_clade=paths[target] is not None,
                         enclosing_hog_parent=paths[target] or '',immediate_gene_tree_parent=parent or '',
                         immediate_parent_descendant_genes=counts[parent] if parent else 1,
                         first_ancestor_with_hog_descendants=join or '',
                         edges_to_first_hog_join=distance if join else '',
                         hogs_below_first_join=hog_counts[join] if join else 0))
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());plan_sha=sha(args.plan)
    def verify():
        if sha(args.plan)!=plan_sha:raise ValueError('Changed plan')
        for path,digest in plan['pins'].items():
            if sha(path)!=digest:raise ValueError('Changed input')
    verify();out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    targets={};flags={};identities={}
    with Path(plan['missing']).open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            label=row['taxon_id']+'_'+row['protein_id'];family=row['family']
            target=flags if row['native_flagged_misplaced']=='True' else targets
            target.setdefault(family,set()).add(label);identities[(family,label)]=row
    csv.field_size_limit(16*1024*1024);parents={f:{} for f in targets}
    with Path(plan['root_hogs']).open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            family=row['OG']
            if family not in targets:continue
            parent=row['Gene Tree Parent Clade']
            if parent in parents[family] or parent=='-':raise ValueError('Ambiguous HOG parent')
            members=set()
            for taxon,cell in row.items():
                if taxon in ('HOG','OG','Gene Tree Parent Clade') or not cell:continue
                members.update(taxon+'_'+protein for protein in cell.split(', '))
            parents[family][parent]=members
    rows=[];seen=set();checked_parents=0
    with Path(plan['resolved_trees']).open() as handle:
        for line in handle:
            family,newick=line.rstrip('\n').split(': ',1)
            if family not in targets:continue
            if family in seen:raise ValueError('Repeated target family tree')
            seen.add(family);checked_parents+=len(parents[family])
            for row in trace_family(newick,parents[family],flags.get(family,set()),targets[family]):
                identity=identities[(family,row['gene_label'])]
                rows.append(dict(family=family,taxon_id=identity['taxon_id'],protein_id=identity['protein_id'],native_gene_id=identity['native_gene_id'],**row))
    if seen!=set(targets) or len(rows)!=sum(map(len,targets.values())):raise ValueError('Incomplete target scope')
    verify();out.mkdir(parents=True);table=out/'gene_tree_positions.tsv'
    with table.open('w') as handle:
        w=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    receipt=dict(status='complete_root_omission_parent_clade_trace',families=len(seen),genes=len(rows),checked_root_hog_parent_memberships=checked_parents,
                 genes_inside_emitted_parent_clades=sum(r['inside_emitted_parent_clade'] for r in rows),
                 genes_outside_emitted_parent_clades=sum(not r['inside_emitted_parent_clade'] for r in rows),
                 genes_with_no_hog_descendant_join=sum(not r['first_ancestor_with_hog_descendants'] for r in rows),
                 plan_sha256=plan_sha,artifacts={'gene_tree_positions.tsv':sha(table)},
                 scope='Every unflagged root omission traced in its resolved tree; all emitted root HOG parent memberships in affected families compared with descendant genes minus native flags. Describes the emitted tree partition, not an independent inference of ancestral homology, duplication timing, loss or contamination.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)


if __name__=='__main__':main()
