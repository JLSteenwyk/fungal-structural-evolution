#!/usr/bin/env python3
"""Measure family, taxon and exact model reuse within all whole-protein candidate screens."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f,delimiter='\t'))


def components(keys, memberships):
    parent = {k:k for k in keys}
    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]; k = parent[k]
        return k
    for members in memberships.values():
        if not members: continue
        anchor = members[0]
        for member in members[1:]:
            parent[find(member)] = find(anchor)
    groups = defaultdict(set)
    for k in keys: groups[find(k)].add(k)
    expected = {frozenset(v) for v in groups.values()}
    # Independent graph traversal on the full incidence graph.
    adjacency = defaultdict(set)
    for members in memberships.values():
        for k in members[1:]:
            adjacency[members[0]].add(k); adjacency[k].add(members[0])
    unseen = set(keys); actual = set()
    while unseen:
        start = unseen.pop(); reached={start}; stack=[start]
        while stack:
            for k in adjacency[stack.pop()]:
                if k in unseen:
                    unseen.remove(k); reached.add(k); stack.append(k)
        actual.add(frozenset(reached))
    assert actual == expected
    return sorted(actual,key=lambda group: sorted(group)[0])


def main():
    cross = Path('results/structural_comparisons/whole-protein-cross-guide-sensitivity-20260927-v1')
    mapping = Path('results/structural_comparisons/whole-protein-common-residues-20260927-v1')
    sources={}
    for root in [cross,mapping]:
        rp=root/'receipt.json';receipt=json.loads(rp.read_text());sources[str(rp)]=sha(rp)
        for name,digest in receipt['artifacts'].items():
            assert sha(root/name)==digest
    triads={r['triad_id']:r for r in rows(mapping/'model_triads.tsv')}
    models=defaultdict(set); link_counts=Counter()
    for r in rows(mapping/'event_reference_triads.tsv'):
        key=(r['family'],r['gene_a'],r['gene_b']); t=triads[r['triad_id']]
        for role in ['a','b','reference']: models[key].add((t[role+'_model'],t[role+'_version']))
        link_counts[key]+=1
    data=rows(cross/'all_pair_comparisons.tsv'); summaries=[]; assignments=[]; frequencies=[]
    for screen in sorted({r['screen'] for r in data}):
        for cohort,flag in [('all_references_eligible','both_guides_all_references_eligible'),('direction_stable','same_structural_direction_both_guides')]:
            selected=[r for r in data if r['screen']==screen and r[flag]=='1']
            keys=[(r['family'],r['gene_a'],r['gene_b']) for r in selected]
            assert len(keys)==len(set(keys))
            families=Counter(k[0] for k in keys);taxa=Counter(k[1].split('_')[0] for k in keys)
            assert all(k[1].split('_')[0]==k[2].split('_')[0] for k in keys)
            incidence=defaultdict(list)
            for k in keys:
                assert models[k]
                for model in models[k]:incidence[('model',)+model].append(k)
            model_components=components(keys,incidence)
            for k in keys:incidence[('family',k[0])].append(k)
            joint_components=components(keys,incidence)
            for kind,groups in [('shared_model',model_components),('family_or_shared_model',joint_components)]:
                for group in groups:
                    cid=hashlib.sha256(json.dumps(sorted(group),separators=(',',':')).encode()).hexdigest()
                    for k in sorted(group):
                        assignments.append(dict(screen=screen,cohort=cohort,dependence=kind,component_id=cid,component_pairs=len(group),family=k[0],gene_a=k[1],gene_b=k[2]))
            model_degrees=[len(v) for k,v in incidence.items() if k[0]=='model']
            summaries.append(dict(screen=screen,cohort=cohort,pairs=len(keys),families=len(families),taxa=len(taxa),
                unique_models=len(model_degrees),models_reused_across_pairs=sum(n>1 for n in model_degrees),
                shared_model_components=len(model_components),largest_shared_model_component=max(map(len,model_components),default=0),
                family_or_model_components=len(joint_components),largest_family_or_model_component=max(map(len,joint_components),default=0),
                largest_family_pairs=max(families.values(),default=0),largest_taxon_pairs=max(taxa.values(),default=0)))
            for axis,counts in [('family',families),('taxon',taxa)]:
                for label,count in sorted(counts.items()): frequencies.append(dict(screen=screen,cohort=cohort,axis=axis,label=label,pairs=count))
    out=Path('results/structural_comparisons/whole-protein-candidate-dependence-20260927-v1');out.mkdir(exist_ok=False);artifacts={}
    for name,data in [('summary.tsv',summaries),('component_membership.tsv',assignments),('sampling_counts.tsv',frequencies)]:
        path=out/name
        with path.open('w') as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');writer.writeheader();writer.writerows(data)
        assert rows(path)==[{k:str(v) for k,v in row.items()} for row in data]
        artifacts[name]=sha(path)
    for path,digest in sources.items():assert sha(path)==digest
    receipt=dict(status='complete_whole_protein_candidate_dependence_inventory',source_hashes=sources,script_sha256=sha(__file__),
        summaries=summaries,artifacts=artifacts,scope='All six screens, eligible and direction-stable subsets. Union of every linked reference and both guides; models identified by accession/version. Union-find and full graph traversal agree for every component. Components identify known reuse, not independent observations: taxa, phylogeny and prediction uncertainty still induce dependence.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
