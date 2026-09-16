#!/usr/bin/env python3
"""Verify Pfam detection states and minimum scores using set-based Fitch traversal."""
import argparse,json
from pathlib import Path
from collections import defaultdict,Counter
from Bio import Phylo
from audit_joint_path_uncertainty import checked,rows,sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['mapping','annotations','domain-inputs','trees','output']:ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=checked(a.mapping);checked(a.annotations);checked(a.domain_inputs)
    for path,digest in r['source_pins'].items():
        if sha(Path(path))!=digest:raise ValueError('Changed pinned source')
    links={(x['marker'],x['taxon_id']):x['sequence_id'] for x in rows(a.domain_inputs/'protein_links.tsv')}
    hits=defaultdict(lambda:defaultdict(list))
    for x in rows(a.annotations/'raw_annotated_hits.tsv'):
        if x['pfam_type'] in ['Domain','Family','Repeat']:hits[x['sequence_id']][x['pfam_accession']].append(float(x['hmm_coverage']))
    overlap={x['sequence_id'] for x in rows(a.annotations/'overlapping_hits.tsv')}
    tips=defaultdict(dict)
    for row in rows(a.mapping/'tip_states.tsv'):
        key=(row['marker'],row['pfam_accession'],row['policy'])
        if row['taxon_id'] in tips[key]:raise ValueError('Duplicate tip state')
        tips[key][row['taxon_id']]=row
    edges=defaultdict(list)
    for row in rows(a.mapping/'possible_change_edges.tsv'):edges[row['marker'],row['pfam_accession'],row['policy']].append(row)
    summaries=list(rows(a.mapping/'profile_summary.tsv'));trees={};scores=tipcount=0;keys=set()
    for row in summaries:
        marker,profile,policy=[row[x] for x in ['marker','pfam_accession','policy']];key=(marker,profile,policy)
        if key in keys:raise ValueError('Duplicate summary')
        keys.add(key)
        if marker not in trees:trees[marker]=Phylo.read(a.trees/marker/'tree.treefile','newick')
        tree=trees[marker];taxa={n.name for n in tree.get_terminals()};states={};counts=Counter()
        for t in taxa:
            sid=links[marker,t];h=hits[sid].get(profile,[])
            if policy=='overlap_partial_unknown' and (sid in overlap or (h and max(h)<.7)):state={0,1};name='unknown'
            elif h:state={1};name='detected'
            else:state={0};name='undetected'
            states[t]=state;counts[name]+=1
        assert all(int(row[k])==counts[k] for k in ['detected','undetected','unknown'])
        eligible=counts['detected']>=5 and counts['undetected']>=5
        assert (row['eligible']=='True')==eligible and int(row['tree_taxa'])==len(taxa)
        if not eligible:
            assert key not in tips and key not in edges and row['minimum_detection_changes']=='';continue
        assert set(tips[key])==taxa
        for t in taxa:
            assert tips[key][t]['sequence_id']==links[marker,t]
            assert {int(c) for c in tips[key][t]['allowed_states']}==states[t];tipcount+=1
        score=0;sets={}
        for node in tree.find_clades(order='postorder'):
            if node.is_terminal():sets[node]=states[node.name];continue
            parts=[sets[c] for c in node.clades];value=parts[0]
            assert len(parts) in [2,3] and (len(parts)==2 or node is tree.root)
            for part in parts[1:]:
                intersection=value&part
                if intersection:value=intersection
                else:value=value|part;score+=1
            sets[node]=value
        assert score==int(row['minimum_detection_changes']);scores+=1
        assert len(edges[key])==int(row['edges_with_change_in_some_optimum'])
        assert sum(x['change_in_every_optimum']=='True' for x in edges[key])==int(row['edges_with_change_in_every_optimum'])
        assert len({x['split_taxa'] for x in edges[key]})==len(edges[key])
    expected={(m,p,policy) for m,tree in trees.items() for p in {p for t in tree.get_terminals() for p in hits[links[m,t.name]]} for policy in ['all_GA_hits','overlap_partial_unknown']}
    assert keys==expected
    result={'status':'passed_marker_pfam_state_and_score_readback','profile_policy_dispositions_checked':len(summaries),'tip_states_checked':tipcount,'independent_fitch_scores_checked':scores,'mapping_receipt_sha256':sha(a.mapping/'receipt.json'),'script_sha256':sha(Path(__file__)),'scope':'All profile dispositions and eligible tip states rebuilt from raw GA hits and overlap table; all minimum scores recomputed with set-based Fitch traversal. Empirical possible/mandatory edge marginals checked for output counts/uniqueness only; numerical correctness of marginals tested separately by exhaustive small-tree enumeration.'}
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
