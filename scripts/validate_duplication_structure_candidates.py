#!/usr/bin/env python3
"""Reconstruct all simple modeled duplication candidates and inspect reported tree nodes."""
import argparse,csv,json,sqlite3
from collections import Counter,defaultdict
from io import StringIO
from pathlib import Path
from Bio import Phylo
from run_ortholog_pair_guide_comparison import sha


def index_tree(newick):
    tree=Phylo.read(StringIO(newick),'newick');nodes={};stack=[tree.root]
    while stack:
        n=stack.pop()
        if not n.name or n.name in nodes:raise ValueError('Unnamed or repeated tree node')
        nodes[n.name]=n;stack.extend(n.clades)
    return nodes


def inspect_node(nodes,name,expected):
    if name not in nodes:return 'missing_reported_node',False
    node=nodes[name];stack=[node];tips=[]
    while stack:
        n=stack.pop()
        if n.clades:stack.extend(n.clades)
        else:tips.append(n.name)
    if len(tips)!=2 or set(tips)!=set(expected):return 'descendant_pair_mismatch',False
    return 'exact_reported_pair',len(node.clades)==2 and all(not n.clades for n in node.clades)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());out=Path(plan['output'])
    def verify():
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed pin '+path)
    verify()
    if out.exists():raise FileExistsError(out)
    source=Path(plan['coverage']);receipt=json.loads((source/'receipt.json').read_text());assert receipt['status']=='complete_duplication_frozen_structure_coverage'
    db=sqlite3.connect('file:'+str(Path(plan['bridge']).resolve())+'?mode=ro',uri=True)
    models={t+'_'+p:(s,m) for t,p,s,m in db.execute('SELECT taxon_id,protein_id,sequence_sha256,model_id FROM structures')};db.close()
    out.mkdir(parents=True);summaries=[];pairs={}
    for entry in plan['guides']:
        tr=json.loads(Path(entry['tree_readback']).read_text())
        if tr['status']!='passed_complete_resolved_tree_membership_readback' or sha(entry['trees'])!=tr['resolved_tree_file_sha256']:raise ValueError('Resolved tree audit mismatch')
        guide=entry['guide'];path=source/(guide+'_terminal_two_model_candidates.tsv')
        assert sha(path)==receipt['artifacts'][path.name]
        with path.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        bykey={(r['family'],r['gene_node']):r for r in rows}
        if len(bykey)!=len(rows):raise ValueError('Repeated candidate event')
        observed=set()
        # Enumerate eligibility anew from every native row and the frozen model lookup.
        with Path(entry['events']).open() as f:
            for event in csv.DictReader(f,delimiter='\t'):
                left,right=event['Genes 1'],event['Genes 2']
                if event['Type']!='Terminal' or not left or not right or ', ' in left or ', ' in right:continue
                if left not in models or right not in models:continue
                ga,gb=sorted([left,right]);sa,ma=models[ga];sb,mb=models[gb]
                key=(event['Orthogroup'],event['Gene Tree Node']);r=bykey.get(key)
                expected=dict(family=key[0],taxon_id=event['Species Tree Node'],gene_node=key[1],support=event['Support'],gene_a=ga,gene_b=gb,model_a=ma,model_b=mb,same_sequence=str(int(sa==sb)),same_model=str(int(ma==mb)))
                if r!=expected or key in observed:raise ValueError('Candidate differs from independent event/model join')
                observed.add(key)
        if observed!=set(bykey):raise ValueError('Missing/extra eligible candidate')
        byfamily=defaultdict(list)
        for r in rows:byfamily[r['family']].append(r)
        checked=set();counts=Counter();pairmap={}
        with (out/(guide+'_candidate_tree_checks.tsv')).open('w') as dst,Path(entry['trees']).open() as src:
            columns=list(rows[0]) + ['tree_status','two_direct_tip_children'];writer=csv.DictWriter(dst,fieldnames=columns,delimiter='\t',lineterminator='\n');writer.writeheader()
            for line in src:
                family,newick=line.rstrip().split(': ',1)
                if family not in byfamily:continue
                if family in checked:raise ValueError('Repeated tree family')
                nodes=index_tree(newick);checked.add(family)
                for r in byfamily[family]:
                    status,direct=inspect_node(nodes,r['gene_node'],[r['gene_a'],r['gene_b']]);writer.writerow(dict(r,tree_status=status,two_direct_tip_children=int(direct)));counts[status]+=1;counts['two_direct_tip_children']+=direct
                    key=(r['gene_a'],r['gene_b'])
                    if key in pairmap:raise ValueError('Repeated candidate protein pair')
                    pairmap[key]=(status,r)
                if len(checked)%1000==0:print(guide,'tree families checked',len(checked),flush=True)
            for family in sorted(set(byfamily)-checked):
                for r in byfamily[family]:
                    writer.writerow(dict(r,tree_status='missing_family_tree',two_direct_tip_children=0));counts['missing_family_tree']+=1
                    key=(r['gene_a'],r['gene_b'])
                    if key in pairmap:raise ValueError('Repeated candidate pair')
                    pairmap[key]=('missing_family_tree',r)
        pairs[guide]=pairmap;summaries.append(dict(guide=guide,candidates=len(rows),independent_candidate_join='passed_all_eligible_native_events',tree_families_checked=len(checked),counts=dict(counts)))
        print(json.dumps(summaries[-1]),flush=True)
    shared=set(pairs['profile'])&set(pairs['mafft']);union=set(pairs['profile'])|set(pairs['mafft'])
    with (out/'guide_pair_membership.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['gene_a','gene_b','profile_tree_status','mafft_tree_status'])
        for key in sorted(union):w.writerow([*key,pairs['profile'].get(key,('not_candidate',))[0],pairs['mafft'].get(key,('not_candidate',))[0]])
    overlap=dict(shared_pairs=len(shared),profile_only=len(set(pairs['profile'])-shared),mafft_only=len(set(pairs['mafft'])-shared),shared_exact_tree_pairs=sum(all(pairs[g][k][0]=='exact_reported_pair' for g in pairs) for k in shared))
    verify();r=dict(status='complete_duplication_candidate_join_and_tree_review',plan_sha256=sha(a.plan),script_sha256=sha(__file__),guides=summaries,overlap=overlap,artifacts={p.name:sha(p) for p in out.iterdir()},scope='Every eligible terminal singleton-side candidate independently rebuilt from native events and frozen models. Reported node descendant pairs checked, with discrepancies retained. Does not independently validate full event-coverage inventory, duplication inference/support, structure quality, rooting, timing or asymmetry; noncandidate and larger events remain outside this candidate review.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n')

if __name__=='__main__':
    csv.field_size_limit(32*1024*1024);main()
