#!/usr/bin/env python3
"""Describe shared sequences, entries and publication links among experimental controls."""
import csv,hashlib,json
from collections import defaultdict,deque
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha

BASE=Path('results/experimental_structures')
KEY=['family','gene_a','gene_b','pfam_accession','screen']


def main():
    sources={}
    def checked(directory,name):
        root=BASE/directory;rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/name
        assert sha(p)==r['artifacts'][name];sources[str(rp)]=sha(rp);sources[str(p)]=sha(p)
        with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
    subjects=checked('whole-domain-case-subject-readback-20260927-v1','subject_metadata_checks.tsv')
    entries={r['entry_id']:r for r in checked('whole-domain-case-metadata-review-20260927-v1','entries.tsv')}
    candidates=checked('whole-domain-case-metadata-links-20260927-v1','annotated_candidate_entities.tsv')
    prior=checked('whole-domain-case-metadata-links-20260927-v1','case_metadata_summary.tsv')
    seq={}
    for r in subjects:
        entity=r['entity_id'];h=r['canonical_sequence_sha256']
        if entity in seq:assert seq[entity]==h
        seq[entity]=h
    selected=sorted(e for e in seq if entries[e.rsplit('_',1)[0]]['methodology']=='experimental')
    parents={e:e for e in selected};buckets=defaultdict(list)
    def find(e):
        while parents[e]!=e:
            parents[e]=parents[parents[e]];e=parents[e]
        return e
    def join(a,b):
        a,b=find(a),find(b)
        if a!=b:parents[b]=a
    for e in selected:
        entry=e.rsplit('_',1)[0];r=entries[entry]
        buckets['sequence',seq[e]].append(e);buckets['entry',entry].append(e)
        doi=r['citation_doi'].strip().lower();pubmed=r['primary_citation_pubmed_id'].strip()
        if doi and doi not in ('unknown','none','?','.'):buckets['doi',doi].append(e)
        if pubmed and pubmed not in ('unknown','None','?','.'):buckets['pubmed',pubmed].append(e)
    edges=[];graph={e:set() for e in selected}
    for (kind,value),members in sorted(buckets.items()):
        for e in members[1:]:
            join(members[0],e);graph[members[0]].add(e);graph[e].add(members[0])
            edges.append(dict(link_type=kind,link_value=value,entity_a=members[0],entity_b=e))
    components=defaultdict(set)
    for e in selected:components[find(e)].add(e)
    # Independent graph traversal checks the full union-find partition.
    unseen=set(selected);traversed=[]
    while unseen:
        start=min(unseen);queue=deque([start]);visited={start}
        while queue:
            for neighbor in graph[queue.popleft()]-visited:visited.add(neighbor);queue.append(neighbor)
        unseen-=visited;traversed.append(frozenset(visited))
    assert set(traversed)=={frozenset(v) for v in components.values()}
    labels={e:hashlib.sha256('\n'.join(sorted(v)).encode()).hexdigest() for v in components.values() for e in v}
    entity_rows=[dict(entity_id=e,entry_id=e.rsplit('_',1)[0],sequence_sha256=seq[e],component_id=labels[e]) for e in selected]
    grouped=defaultdict(list)
    for r in candidates:grouped[tuple(r[k] for k in KEY)].append(r)
    summary=[]
    for r in prior:
        group=grouped[tuple(r[k] for k in KEY)];field='shared_domain_pass' if r['coverage_region']=='domain' else 'shared_domain_and_outside_pass'
        entities={x['entity_id'] for x in group if x[field]=='1' and x['methodology']=='experimental'}
        assert len(entities)==int(r['experimental_entities'])
        sizes=defaultdict(int)
        for e in entities:sizes[labels[e]]+=1
        summary.append(dict(r,distinct_exact_sequences=len({seq[e] for e in entities}),linked_components=len(sizes),largest_component_entities=max(sizes.values(),default=0)))
    out=BASE/'whole-domain-case-experimental-dependence-20260927-v1';out.mkdir(exist_ok=False);artifacts={}
    for name,data in [('entity_components.tsv',entity_rows),('dependency_edges.tsv',edges),('case_summary.tsv',summary)]:
        p=out/name
        with p.open('w') as f:w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
        with p.open() as f:assert list(csv.DictReader(f,delimiter='\t'))==[{k:str(v) for k,v in row.items()} for row in data]
        artifacts[name]=sha(p)
    for p,digest in sources.items():assert sha(p)==digest
    result=dict(status='complete_case_experimental_dependence_census',source_hashes=sources,script_sha256=sha(__file__),experimental_entities=len(selected),exact_sequences=len({seq[e] for e in selected}),components=len(components),largest_component=max(map(len,components.values())),case_screen_region_rows=len(summary),artifacts=artifacts,scope='Conservative linkage census on all experimentally classified candidate entities, before coverage selection. Edges represent exact canonical sequence, shared entry, DOI or PubMed identifier; unknown citations never joined. Full union-find components independently matched graph traversal. Components can connect different proteins through a paper or entry and are not proven independent experiments, phylogenetic units or effective sample sizes. No structures discarded; no-hits remain in case summaries.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))
    for r in summary:
        if r['screen']=='n50_c70' and r['coverage_region']=='domain_and_outside':print(r['family'],r['experimental_entities'],r['distinct_exact_sequences'],r['linked_components'],r['largest_component_entities'])


if __name__=='__main__':main()
