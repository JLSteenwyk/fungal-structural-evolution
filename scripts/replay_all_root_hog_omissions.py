#!/usr/bin/env python3
"""Replay every affected resolved family with original native gene identifiers."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import time
from orthofinder.gene_tree_inference import trees2ologs_of as native
from replay_root_hog_omission_case import sha


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());plan_hash=sha(args.plan)
    def verify():
        if sha(args.plan)!=plan_hash:raise ValueError('Plan changed')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Input changed: '+p)
    verify();out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    targets=defaultdict(set);flags=defaultdict(set);sizes={}
    with Path(plan['missing']).open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            label=r['taxon_id']+'_'+r['protein_id'];fam=r['family'];sizes[fam]=int(r['family_size'])
            (flags if r['native_flagged_misplaced']=='True' else targets)[fam].add(label)
    species={}
    for line in Path(plan['species_ids']).read_text().splitlines():
        sid,name=line.split(': ',1);species[sid]=name.rsplit('.',1)[0]
    reverse={v:k for k,v in species.items()}
    st=native.tree_lib.Tree(plan['species_tree'],format=1)
    for leaf in st:leaf.name=reverse[leaf.name]
    neighbours=native.GetSpeciesNeighbours(st)
    identities={}
    with Path(plan['sequence_ids']).open() as f:
        for line in f:
            gid,protein=line.rstrip('\n').split(': ',1)
            identities[species[gid.split('_')[0]]+'_'+protein]=gid
    csv.field_size_limit(16*1024*1024);expected={fam:{} for fam in targets}
    with Path(plan['root_hogs']).open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            fam=row['OG']
            if fam not in targets:continue
            parent=row['Gene Tree Parent Clade']
            if parent in expected[fam]:raise ValueError('Duplicate parent')
            expected[fam][parent]={taxon+'_'+p for taxon in species.values() for p in row[taxon].split(', ') if p}
    trees={}
    with Path(plan['resolved_trees']).open() as f:
        for line in f:
            fam,newick=line.rstrip('\n').split(': ',1)
            if fam in targets:
                if fam in trees:raise ValueError('Duplicate tree')
                trees[fam]=newick
    if set(trees)!=set(targets):raise ValueError('Incomplete affected trees')
    counts=defaultdict(int);started=time.monotonic()
    with (out/'family_replays.jsonl').open('x') as log:
        for fam in sorted(targets,key=lambda f:(sizes[f],f)):
            row={'family':fam,'source_genes':sizes[fam],'unflagged_missing_genes':len(targets[fam])}
            try:
                tree=native.tree_lib.Tree(trees.pop(fam),format=1);labels={}
                for leaf in tree:
                    original=leaf.name;leaf.name=identities[original];labels[leaf.name]=original
                if len(labels)!=sizes[fam]:raise ValueError('Source/tree size mismatch')
                _,tree,suspect,dups=native.GetOrthologues_from_tree(int(fam[2:]),tree,st,native.GeneToSpecies_dash,neighbours,q_get_dups=True,qNoRecon=True)
                writer=native.HogWriter(st,[n.name for n in st.traverse() if not n.is_leaf()],labels,species,list(map(int,species)),False,write_output=False)
                writer.mark_dups_below(tree);replay={};states={}
                for node in tree.traverse('preorder'):
                    if not node.is_leaf():states[node.name]={'sp_node':node.sp_node,'dup':node.dup,'dups_below':sorted(node.dups_below)}
                    for level,r in writer.write_clade_v2(node,fam):
                        if level=='N0':
                            if r[1] in replay:raise ValueError('Repeated replay parent')
                            replay[r[1]]={g for cell in r[2:] for g in cell.split(', ') if g}
                replay_flags={labels[g] for g in suspect};omitted=set(labels.values())-set().union(set(),*replay.values())
                root_equal=replay==expected[fam];flag_equal=replay_flags==flags.get(fam,set())
                missing_equal=omitted==targets[fam]|flags.get(fam,set())
                row.update(status='matched' if root_equal and flag_equal and missing_equal else 'mismatch',root_memberships_match=root_equal,misplaced_flags_match=flag_equal,omissions_match=missing_equal,root_groups=len(replay),replayed_flagged_genes=len(replay_flags),replayed_missing_genes=len(omitted))
                positions=[]
                for leaf in tree:
                    if labels[leaf.name] in targets[fam]:
                        parent=leaf.up
                        positions.append({'gene':labels[leaf.name],'parent':parent.name,**states[parent.name]})
                row['omitted_gene_parent_states']=positions
                if row['status']=='mismatch':
                    row['differing_parent_clades']=sorted(k for k in set(replay)|set(expected[fam]) if replay.get(k)!=expected[fam].get(k))
                    row['flag_difference']=sorted(replay_flags.symmetric_difference(flags.get(fam,set())))
                    row['omission_difference']=sorted(omitted.symmetric_difference(targets[fam]|flags.get(fam,set())))
            except Exception as e:
                row.update(status='error',error_type=type(e).__name__,error=str(e))
            counts[row['status']]+=1;counts['processed']+=1
            log.write(json.dumps(row,sort_keys=True)+'\n');log.flush()
            (out/'state.json').write_text(json.dumps({'families_total':len(targets),'counts':dict(counts),'current_family':fam,'elapsed_seconds':time.monotonic()-started},indent=2)+'\n')
    verify()
    receipt={'status':'passed_all_affected_family_replays' if counts['matched']==len(targets) else 'completed_with_replay_discrepancies','families':len(targets),'unflagged_missing_genes':sum(map(len,targets.values())),'counts':dict(counts),'plan_sha256':plan_hash,'artifacts':{'family_replays.jsonl':sha(out/'family_replays.jsonl')},'scope':'Native classification replay from resolved trees with original gene IDs; exact root parent/membership, flagged identities and omission comparisons. Preserves discrepancies. Not independent rooting, initial resolution, or biological event validation.'}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()
