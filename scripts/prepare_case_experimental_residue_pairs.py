#!/usr/bin/env python3
"""Freeze exact query/entity residue correspondences for all case search contexts."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from screen_duplication_domain_alignment_coverage import sha

BASE=Path('results')
SEARCH=BASE/'experimental_structures/whole-domain-case-sequence-search-20260927-v1'
METADATA=BASE/'experimental_structures/whole-domain-case-metadata-20260927-v1'
AUDIT=BASE/'experimental_structures/whole-domain-case-subject-readback-20260927-v1'
FASTA=BASE/'structural_comparisons/case-independent-control-inputs-20260927-v1/missing_independent_sequences.faa'
OUTPUT=BASE/'experimental_structures/whole-domain-case-residue-pairs-20260927-v1'


def main():
    audit=json.loads((AUDIT/'receipt.json').read_text())
    assert audit['status']=='complete_case_subject_metadata_readback'
    sources={str(AUDIT/'receipt.json'):sha(AUDIT/'receipt.json')}
    for path in [SEARCH/'receipt.json',SEARCH/'config.json',SEARCH/'sequence_dispositions.json',METADATA/'receipt.json',FASTA]:sources[str(path)]=sha(path)
    sr=json.loads((SEARCH/'receipt.json').read_text());config=json.loads((SEARCH/'config.json').read_text())
    assert sr['config_sha256']==sha(SEARCH/'config.json') and config['fasta_sha256']==sha(FASTA)
    assert sr['artifacts']['sequence_dispositions.json']==sha(SEARCH/'sequence_dispositions.json')
    assert audit['source_hashes'][str(METADATA/'receipt.json')]==sha(METADATA/'receipt.json')
    queries={r.id:str(r.seq) for r in SeqIO.parse(FASTA,'fasta')}
    subjects={}
    for r in json.loads((METADATA/'receipt.json').read_text())['responses']:
        if r['kind']!='polymer_entity':continue
        path=METADATA/'polymer_entity'/(r['identifier']+'.json');assert sha(path)==r['response_sha256'];sources[str(path)]=sha(path)
        subjects[r['identifier']]=''.join(json.loads(path.read_text())['entity_poly']['pdbx_seq_one_letter_code_can'].split())
    OUTPUT.mkdir(exist_ok=False);path=OUTPUT/'residue_correspondences.jsonl.gz'
    contexts=paired_total=identity_total=query_gap_total=subject_gap_total=0;keys=set();dispositions=[]
    with gzip.open(path,'wt') as output:
        for d in json.loads((SEARCH/'sequence_dispositions.json').read_text()):
            sid=d['sequence_id'];query=queries[sid];rp=SEARCH/(sid+'.json');assert sha(rp)==d['response_sha256'];sources[str(rp)]=sha(rp)
            hits=json.loads(rp.read_text()).get('result_set',[]);assert len(hits)==d['hits'];sequence_contexts=0
            for hit in hits:
                entity=hit['identifier'];subject=subjects[entity]
                matches=[x for service in hit['services'] for node in service['nodes'] for x in node['match_context']]
                for index,c in enumerate(matches):
                    key=(sid,entity,index);assert key not in keys;keys.add(key)
                    qa,sa=c['query_aligned_seq'],c['subject_aligned_seq']
                    assert qa.replace('-','')==query[c['query_beg']-1:c['query_end']]
                    assert sa.replace('-','')==subject[c['subject_beg']-1:c['subject_end']]
                    qpos,spos=c['query_beg']-1,c['subject_beg']-1
                    pairs=[];identical=[];query_only=[];subject_only=[]
                    for q,s in zip(qa,sa):
                        assert q!='-' or s!='-'
                        if q!='-':qpos+=1
                        if s!='-':spos+=1
                        if q!='-' and s!='-':pairs.append([qpos,spos]);identical.append(int(q==s))
                        elif q!='-':query_only.append(qpos)
                        else:subject_only.append(spos)
                    assert len(qa)==len(sa)==c['alignment_length'] and qpos==c['query_end'] and spos==c['subject_end']
                    # Independent column-to-position lookup verifies offsets and each pair.
                    qi={col:c['query_beg']+i for i,col in enumerate(j for j,x in enumerate(qa) if x!='-')}
                    si={col:c['subject_beg']+i for i,col in enumerate(j for j,x in enumerate(sa) if x!='-')}
                    assert pairs==[[qi[col],si[col]] for col in sorted(qi.keys()&si.keys())]
                    assert query_only==[qi[col] for col in sorted(qi.keys()-si.keys())]
                    assert subject_only==[si[col] for col in sorted(si.keys()-qi.keys())]
                    assert identical==[int(query[q-1]==subject[s-1]) for q,s in pairs]
                    record=dict(sequence_id=sid,entity_id=entity,context_index=index,query_length=len(query),subject_length=len(subject),query_sha256=hashlib.sha256(query.encode()).hexdigest(),subject_sha256=hashlib.sha256(subject.encode()).hexdigest(),query_subject_position_pairs=pairs,identical=identical,query_positions_aligned_to_gap=query_only,subject_positions_aligned_to_gap=subject_only)
                    encoded=json.dumps(record,separators=(',',':'));assert json.loads(encoded)==record;output.write(encoded+'\n')
                    contexts+=1;sequence_contexts+=1;paired_total+=len(pairs);identity_total+=sum(identical);query_gap_total+=len(query_only);subject_gap_total+=len(subject_only)
            dispositions.append(dict(sequence_id=sid,hits=len(hits),contexts=sequence_contexts))
    assert contexts==audit['alignment_contexts']==3091 and len(queries)==len(dispositions)==39
    with gzip.open(path,'rt') as handle:
        actual=[(r['sequence_id'],r['entity_id'],r['context_index']) for r in map(json.loads,handle)]
    assert len(actual)==contexts and set(actual)==keys
    dp=OUTPUT/'sequence_dispositions.json';dp.write_text(json.dumps(dispositions,indent=2)+'\n')
    for p,digest in sources.items():assert sha(p)==digest
    result=dict(status='complete_all_case_experimental_residue_correspondences',source_hashes=sources,script_sha256=sha(__file__),sequences=39,contexts=contexts,paired_residue_occurrences=paired_total,identical_residue_occurrences=identity_total,query_residues_aligned_to_gap=query_gap_total,subject_residues_aligned_to_gap=subject_gap_total,artifacts={p.name:sha(p) for p in [path,dp]},scope='All search contexts retained, including integrative candidates. One-based query/entity coordinates and insertion/deletion positions reconstructed against exact canonical sequences; independently checked by alignment-column indexing. Paired occurrences are not independent observations, observed experimental residues, or a common-residue fit. No structural or functional conclusion.')
    (OUTPUT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
