#!/usr/bin/env python3
"""Check every retrieved entity/entry and every search subject against canonical metadata."""
import csv,hashlib,json,subprocess,time
from collections import Counter
from pathlib import Path
import psutil
from screen_duplication_domain_alignment_coverage import sha


def main():
    lp=Path('metadata/case_experimental_metadata_launch_20260927.json');launch=json.loads(lp.read_text());lh=sha(lp)
    while True:
        try:
            p=psutil.Process(launch['pid'])
            if p.create_time()!=launch['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        print('Waiting for exact case metadata retrieval',launch['pid'],flush=True);time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0') and sha(lp)==lh
    assert sha('scripts/retrieve_experimental_metadata.py')==launch['script_sha256']
    base=Path('results/experimental_structures');root=base/'whole-domain-case-metadata-20260927-v1';rp=root/'receipt.json';r=json.loads(rp.read_text());cp=root/'config.json';c=json.loads(cp.read_text())
    assert r['status']=='complete_experimental_candidate_metadata' and r['config_sha256']==sha(cp)
    assert c['inventory_receipt_sha256']==launch['inventory_receipt_sha256']==sha(base/'whole-domain-case-metadata-inputs-20260927-v1/receipt.json')
    ids=set((base/'whole-domain-case-metadata-inputs-20260927-v1/polymer_entity_ids.txt').read_text().split());entries={i.rsplit('_',1)[0] for i in ids}
    expected={('polymer_entity',i) for i in ids}|{('entry',i) for i in entries};seen=set();data={};hashes={str(p):sha(p) for p in [lp,rp,cp]}
    for row in r['responses']:
        k=row['kind'],row['identifier'];assert k in expected and k not in seen;seen.add(k)
        path=root/k[0]/(k[1]+'.json');receipt=path.with_suffix('.receipt.json');rr=json.loads(receipt.read_text())
        assert sha(path)==row['response_sha256']==rr['response_sha256'] and sha(receipt)==row['receipt_sha256'] and rr['config_sha256']==sha(cp)
        obj=json.loads(path.read_text());assert obj['rcsb_id']==k[1] and rr['identifier']==k[1] and rr['kind']==k[0]
        suffix=k[1].replace('_','/') if k[0]=='polymer_entity' else k[1];assert rr['url']==c['endpoint']+'/'+k[0]+'/'+suffix
        data[k]=obj;hashes[str(path)]=sha(path);hashes[str(receipt)]=sha(receipt)
    assert seen==expected and r['entity_count']==len(ids)==1205 and r['entry_count']==len(entries)==707
    search=base/'whole-domain-case-sequence-search-20260927-v1';sr=search/'receipt.json';s=json.loads(sr.read_text());hashes[str(sr)]=sha(sr)
    dp=search/'sequence_dispositions.json';assert sha(dp)==s['artifacts'][dp.name];hashes[str(dp)]=sha(dp)
    output=[]
    for disposition in json.loads(dp.read_text()):
        sid=disposition['sequence_id'];p=search/(sid+'.json');assert sha(p)==disposition['response_sha256'];hashes[str(p)]=sha(p)
        for hit in json.loads(p.read_text()).get('result_set',[]):
            entity=hit['identifier'];obj=data['polymer_entity',entity];entry=data['entry',entity.rsplit('_',1)[0]]
            seq=''.join(obj['entity_poly']['pdbx_seq_one_letter_code_can'].split());assert seq
            contexts=[ctx for service in hit['services'] for node in service['nodes'] for ctx in node['match_context']]
            for index,ctx in enumerate(contexts):
                lo=ctx['subject_beg'];hi=ctx['subject_end'];subject=ctx['subject_aligned_seq'].replace('-','')
                length_ok=len(seq)==ctx['subject_length'];match=seq[lo-1:hi]==subject
                status='exact_canonical_subject_match' if length_ok and match else 'subject_metadata_mismatch_requires_review'
                output.append(dict(sequence_id=sid,entity_id=entity,context_index=index,status=status,canonical_length=len(seq),search_subject_length=ctx['subject_length'],subject_begin=lo,subject_end=hi,canonical_sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),polymer_type=obj['entity_poly'].get('type',''),reported_mutation_count=obj['entity_poly'].get('rcsb_mutation_count',''),reported_nonstandard_monomers=obj['entity_poly'].get('rcsb_non_std_monomer_count',''),experimental_methods_json=json.dumps([x['method'] for x in entry.get('exptl',[])],separators=(',',':')),resolution_combined_json=json.dumps(entry.get('rcsb_entry_info',{}).get('resolution_combined',[]),separators=(',',':')),release_date=entry.get('rcsb_accession_info',{}).get('initial_release_date','')))
    assert len(output)==s['sequence_entity_hits']==3091
    out=base/'whole-domain-case-subject-readback-20260927-v1';out.mkdir(exist_ok=False);p=out/'subject_metadata_checks.tsv'
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]),delimiter='\t');w.writeheader();w.writerows(output)
    with p.open() as f:assert list(csv.DictReader(f,delimiter='\t'))==[{k:str(v) for k,v in row.items()} for row in output]
    for path,digest in hashes.items():assert sha(path)==digest
    result=dict(status='complete_case_subject_metadata_readback',terminal_state=state,source_hashes=hashes,checker_sha256=sha(__file__),entities=len(ids),entries=len(entries),alignment_contexts=len(output),dispositions=dict(Counter(x['status'] for x in output)),artifacts={p.name:sha(p)},scope='All metadata response identities/hashes and subject alignment substrings checked. Mismatches explicit; mutation, method, resolution and release annotations exported without selecting favorable entries. Canonical entity sequence is not observed residue coverage; coordinate and construct validation remain outstanding.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
