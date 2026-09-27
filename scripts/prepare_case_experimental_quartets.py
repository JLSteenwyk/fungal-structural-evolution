#!/usr/bin/env python3
"""Intersect all three fungal/experimental sequence mappings without selecting hits."""
import csv
import gzip
import itertools
import json
from collections import defaultdict
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha

BASE=Path('results')
OUTPUT=BASE/'experimental_structures/whole-domain-case-sequence-quartets-20260927-v1'
ROLES=['a','b','reference']


def table(path):
    with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def main():
    sources={}
    def checked(root,name):
        rp=root/'receipt.json';r=json.loads(rp.read_text());path=root/name
        assert sha(path)==r['artifacts'][name]
        sources[str(rp)]=sha(rp);sources[str(path)]=sha(path)
        return path
    structural=BASE/'structural_comparisons'
    links=table(checked(structural/'whole-domain-case-dossiers-20260927-v1','whole_reference_links.tsv'))
    wanted={r['triad_id'] for r in links}
    triads={r['triad_id']:r for r in table(checked(structural/'whole-protein-common-residues-20260927-v1','model_triads.tsv')) if r['triad_id'] in wanted}
    assert set(triads)==wanted and len(triads)==13
    root=BASE/'experimental_structures/whole-domain-case-residue-pairs-20260927-v1'
    mappings=defaultdict(lambda:defaultdict(list))
    with gzip.open(checked(root,'residue_correspondences.jsonl.gz'),'rt') as handle:
        for r in map(json.loads,handle):mappings[r['sequence_id']][r['entity_id']].append(r)
    dispositions=json.loads(checked(root,'sequence_dispositions.json').read_text())
    available={r['sequence_id'] for r in dispositions}
    OUTPUT.mkdir(exist_ok=False);path=OUTPUT/'quartet_residue_maps.jsonl.gz'
    counts=[];n=positions_total=0;expected_keys=set()
    with gzip.open(path,'wt') as output:
        for tid,triad in sorted(triads.items()):
            sequences=['S'+triad[role+'_sequence_sha256'] for role in ROLES]
            assert set(sequences)<=available
            entity_sets=[set(mappings[s]) for s in sequences]
            common=set.intersection(*entity_sets);union=set.union(*entity_sets)
            count=0;empty=0
            for entity in sorted(common):
                for contexts in itertools.product(*(mappings[s][entity] for s in sequences)):
                    inverse=[{s:q for q,s in r['query_subject_position_pairs']} for r in contexts]
                    assert all(len(d)==len(r['query_subject_position_pairs']) for d,r in zip(inverse,contexts))
                    shared=sorted(set.intersection(*(set(d) for d in inverse)))
                    quartets=[[d[s] for d in inverse]+[s] for s in shared]
                    # Independent ordered membership scan must produce the exact intersection.
                    alternate=[]
                    lookup_b=dict((s,q) for q,s in contexts[1]['query_subject_position_pairs'])
                    lookup_r=dict((s,q) for q,s in contexts[2]['query_subject_position_pairs'])
                    for q,s in contexts[0]['query_subject_position_pairs']:
                        if s in lookup_b and s in lookup_r:alternate.append([q,lookup_b[s],lookup_r[s],s])
                    assert alternate==quartets
                    assert all(len({row[i] for row in quartets})==len(quartets) for i in range(4))
                    assert len({r['subject_sha256'] for r in contexts})==1
                    indices=[r['context_index'] for r in contexts]
                    key=(tid,entity,tuple(indices));assert key not in expected_keys;expected_keys.add(key)
                    record=dict(triad_id=tid,entity_id=entity,context_indices=indices,
                                sequence_ids=sequences,model_ids=[triad[r+'_model'] for r in ROLES],
                                model_versions=[int(triad[r+'_version']) for r in ROLES],
                                query_lengths=[r['query_length'] for r in contexts],
                                subject_length=contexts[0]['subject_length'],subject_sha256=contexts[0]['subject_sha256'],
                                position_order=['a','b','reference','experimental_entity'],
                                common_positions=quartets,common_count=len(quartets))
                    encoded=json.dumps(record,separators=(',',':'));assert json.loads(encoded)==record
                    output.write(encoded+'\n');count+=1;n+=1;positions_total+=len(quartets);empty+=not quartets
            counts.append(dict(triad_id=tid,sequence_ids=sequences,role_entity_counts=[len(x) for x in entity_sets],union_entity_count=len(union),shared_entity_count=len(common),context_combinations=count,empty_common_maps=empty,disposition='shared_sequence_candidates' if common else 'no_shared_sequence_candidate'))
    with gzip.open(path,'rt') as handle:
        records=list(map(json.loads,handle))
    assert len(records)==n and {(r['triad_id'],r['entity_id'],tuple(r['context_indices'])) for r in records}==expected_keys
    assert sum(len(r['common_positions']) for r in records)==positions_total
    cp=OUTPUT/'triad_dispositions.json';cp.write_text(json.dumps(counts,indent=2)+'\n')
    lp=OUTPUT/'case_reference_links.tsv'
    with lp.open('w') as handle:
        writer=csv.DictWriter(handle,list(links[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(links)
    assert table(lp)==links
    for p,digest in sources.items():assert sha(p)==digest
    result=dict(status='complete_case_experimental_sequence_quartets',source_hashes=sources,script_sha256=sha(__file__),triads=len(triads),triads_with_shared_candidates=sum(r['shared_entity_count']>0 for r in counts),quartet_maps=n,common_residue_occurrences=positions_total,artifacts={p.name:sha(p) for p in [path,cp,lp]},scope='Common sequence correspondences anchored to one experimental entity for A/B/reference. Every common hit and context combination retained, plus explicit no-shared-hit triads. Sequence-only definition differs from structural alignment maps; coordinate availability, domain partitions, confidence masks, methodology exclusion, fit rank and biological interpretation remain downstream.')
    (OUTPUT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))
    for r in counts:print(r['triad_id'][:10],r['role_entity_counts'],r['shared_entity_count'])


if __name__=='__main__':main()
