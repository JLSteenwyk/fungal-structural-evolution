#!/usr/bin/env python3
"""Map minimum changes in detectable Pfam profiles, retaining annotation uncertainty."""
import argparse,json
from pathlib import Path
from collections import defaultdict,Counter
from Bio import Phylo,SeqIO
from audit_joint_path_uncertainty import checked,rows,sha
from prepare_paired_phylogenetic_inputs import write_table


def sankoff(tree,states):
    nodes=list(tree.find_clades(order='postorder'));down={};outside={tree.root:[0,0]};contribution={}
    for node in nodes:
        if node.is_terminal():down[node]=[0 if i in states[node.name] else 10**9 for i in [0,1]]
        else:
            down[node]=[0,0]
            for child in node.clades:
                contribution[child]=[min(down[child][t]+(s!=t) for t in [0,1]) for s in [0,1]]
                down[node]=[down[node][s]+contribution[child][s] for s in [0,1]]
    score=min(down[tree.root]);edge=[]
    for parent in tree.find_clades(order='preorder'):
        for child in parent.clades:
            rest=[outside[parent][s]+sum(contribution[c][s] for c in parent.clades if c is not child) for s in [0,1]]
            outside[child]=[min(rest[s]+(s!=t) for s in [0,1]) for t in [0,1]]
            pairs=[(s,t) for s in [0,1] for t in [0,1] if rest[s]+(s!=t)+down[child][t]==score]
            if not pairs:raise ValueError('No globally optimal edge state pair')
            could=any(s!=t for s,t in pairs);must=all(s!=t for s,t in pairs)
            edge.append((child,could,must))
    return score,edge


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['annotations','domain-inputs','snapshot','readback','trees','output']:ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    ar=checked(a.annotations);ir=checked(a.domain_inputs);sr=checked(a.snapshot);rr=checked(a.readback)
    if rr['snapshot_receipt_sha256']!=sha(a.snapshot/'receipt.json'):raise ValueError('Snapshot readback mismatch')
    links={};sequences={r.id:str(r.seq) for r in SeqIO.parse(a.domain_inputs/'sequences.faa','fasta')}
    for row in rows(a.domain_inputs/'protein_links.tsv'):
        key=(row['marker'],row['taxon_id'])
        if key in links:raise ValueError('Multiple selected proteins per marker/taxon')
        links[key]=row['sequence_id']
    hits=defaultdict(lambda:defaultdict(list));meta={}
    for row in rows(a.annotations/'raw_annotated_hits.tsv'):
        if row['pfam_type'] in ['Domain','Family','Repeat']:
            hits[row['sequence_id']][row['pfam_accession']].append(row)
            meta[row['pfam_accession']]=(row['pfam_name'],row['pfam_type'])
    ambiguous={x['sequence_id'] for x in rows(a.annotations/'overlapping_hits.tsv')}
    pins={str(p/'receipt.json'):sha(p/'receipt.json') for p in [a.annotations,a.domain_inputs,a.snapshot,a.readback]}
    summaries=[];tiprows=[];edgerows=[]
    for source in sr['inputs']:
        marker=source['marker'];path=a.trees/marker/'tree.treefile'
        if sha(path)!=source['tree_sha256']:raise ValueError('Changed tree')
        pins[str(path)]=sha(path);tree=Phylo.read(path,'newick');taxa={x.name for x in tree.get_terminals()}
        mapping={t:links[marker,t] for t in taxa};profiles=sorted({p for sid in mapping.values() for p in hits[sid]})
        alignment={r.id:str(r.seq) for r in SeqIO.parse(a.trees/marker/'input.faa','fasta')}
        for t,sid in mapping.items():
            # Trimmed marker alignment must preserve order of observed source residues.
            it=iter(sequences[sid]);assert all(c in it for c in alignment[t] if c in 'ACDEFGHIKLMNPQRSTVWY')
        splits={}
        for node in tree.find_clades():
            if node is tree.root:continue
            side={t.name for t in node.get_terminals()};split=min([tuple(sorted(side)),tuple(sorted(taxa-side))],key=lambda x:(len(x),x));splits[node]=';'.join(split)
        for profile in profiles:
            for policy in ['all_GA_hits','overlap_partial_unknown']:
                states={};counts=Counter()
                for t,sid in mapping.items():
                    ph=hits[sid].get(profile,[])
                    unknown=policy=='overlap_partial_unknown' and (sid in ambiguous or (ph and not any(float(h['hmm_coverage'])>=.7 for h in ph)))
                    value='unknown' if unknown else 'detected' if ph else 'undetected'
                    states[t]={0,1} if unknown else {1} if ph else {0};counts[value]+=1
                eligible=counts['detected']>=5 and counts['undetected']>=5
                score='';could=must=0
                if eligible:
                    score,edges=sankoff(tree,states)
                    for t in sorted(taxa):tiprows.append(dict(marker=marker,pfam_accession=profile,policy=policy,taxon_id=t,sequence_id=mapping[t],allowed_states=''.join(map(str,sorted(states[t])))))
                    for node,c,m in edges:
                        if c:
                            could+=1;must+=int(m)
                            edgerows.append(dict(marker=marker,pfam_accession=profile,policy=policy,split_taxa=splits[node],terminal=node.is_terminal(),change_in_every_optimum=m,sh_alrt_percent='' if node.is_terminal() or node.confidence is None else node.confidence))
                summaries.append(dict(marker=marker,pfam_accession=profile,pfam_name=meta[profile][0],pfam_type=meta[profile][1],policy=policy,tree_taxa=len(taxa),detected=counts['detected'],undetected=counts['undetected'],unknown=counts['unknown'],eligible=eligible,minimum_detection_changes=score,edges_with_change_in_some_optimum=could,edges_with_change_in_every_optimum=must))
        print(marker,'complete',flush=True)
    a.output.mkdir(parents=True)
    for name,data in [('profile_summary.tsv',summaries),('tip_states.tsv',tiprows),('possible_change_edges.tsv',edgerows)]:write_table(a.output/name,data)
    result=dict(status='complete_marker_pfam_detection_parsimony',markers=sr['completed_markers'],planned_markers=sr['planned_markers'],profile_policy_rows=len(summaries),eligible_profile_policy_rows=sum(x['eligible'] for x in summaries),possible_change_edge_rows=len(edgerows),source_pins=pins,script_sha256=sha(Path(__file__)),interpretation='Symmetric unit-cost binary Sankoff on unrooted marker trees. Zero is no detectable GA profile, not confirmed domain absence. Strict policy masks all sequence-level overlap cases and partial-only profile hits below 70% HMM coverage. At least five detected and five undetected tips required; all profiles retained in dispositions. No rooted gain/loss polarity, copy-number/fusion/rearrangement inference, calibrated confidence, species reconciliation or causal claim. Every possible/mandatory edge is conditional on this annotation and minimum-change criterion; unsupported trees and missing-marker bias remain.',artifacts={p.name:sha(p) for p in a.output.iterdir()})
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_pins','artifacts']},indent=2))

if __name__=='__main__':main()
