#!/usr/bin/env python3
"""Cache complete experimental RCSB sequence-search responses for every case sequence."""
import datetime,fcntl,hashlib,json,time,urllib.request,urllib.parse
from pathlib import Path
from Bio import SeqIO
from screen_duplication_domain_alignment_coverage import sha


def main():
    source=Path('results/structural_comparisons/case-independent-control-inputs-20260927-v1');sr=source/'receipt.json';s=json.loads(sr.read_text())
    fasta=source/'missing_independent_sequences.faa';assert sha(fasta)==s['artifacts'][fasta.name]
    sequences={r.id:str(r.seq) for r in SeqIO.parse(fasta,'fasta')};assert len(sequences)==39
    assert all(sid=='S'+hashlib.sha256(seq.encode()).hexdigest() for sid,seq in sequences.items())
    out=Path('results/experimental_structures/whole-domain-case-sequence-search-20260927-v1');out.mkdir(exist_ok=True)
    lock=(out/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    config=dict(source_receipt_sha256=sha(sr),fasta_sha256=sha(fasta),script_sha256=sha(__file__),endpoint='https://search.rcsb.org/rcsbsearch/v2/query',identity_cutoff=.3,evalue_cutoff=1e-5,content_type='experimental',results_verbosity='verbose',resources=dict(cpus=1,http_workers=1,memory_gib=2,output_gib=.5,planning_minutes=[2,30],paid_resources=False),scope='Candidate experimental homolog nomination for all 39 sequences; retain all hits without representative selection or coverage filtering. Exact sequence, construct, domain coverage, resolution and training overlap remain downstream.')
    cp=out/'config.json'
    if cp.exists():assert json.loads(cp.read_text())==config
    else:cp.write_text(json.dumps(config,indent=2)+'\n')
    dispositions=[];entities=set()
    for number,(sid,seq) in enumerate(sorted(sequences.items()),1):
        query=dict(query=dict(type='terminal',service='sequence',parameters=dict(sequence_type='protein',value=seq,identity_cutoff=.3,evalue_cutoff=1e-5)),return_type='polymer_entity',request_options=dict(return_all_hits=True,results_content_type=['experimental'],results_verbosity='verbose',scoring_strategy='sequence'))
        rp=out/(sid+'.receipt.json');responsepath=out/(sid+'.json')
        if rp.exists():
            record=json.loads(rp.read_text());assert record['query']==query and record['config_sha256']==sha(cp) and record['response_sha256']==sha(responsepath)
        else:
            url=config['endpoint']+'?'+urllib.parse.urlencode({'json':json.dumps(query)})
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(url,timeout=45) as response:status=response.status;raw=response.read()
                    assert status in [200,204]
                    payload=json.loads(raw) if status==200 else dict(total_count=0,result_set=[])
                    hits=payload.get('result_set',[]);assert payload['total_count']==len(hits)
                    assert len({h['identifier'] for h in hits})==len(hits)
                    break
                except Exception:
                    if attempt==3:raise
                    time.sleep(2**attempt)
            responsepath.write_bytes(raw if status==200 else json.dumps(payload).encode())
            record=dict(sequence_id=sid,query=query,config_sha256=sha(cp),response_sha256=sha(responsepath),http_status=status,total_count=len(hits),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
            rp.write_text(json.dumps(record,indent=2)+'\n');time.sleep(.25)
        payload=json.loads(responsepath.read_text());hits=payload.get('result_set',[])
        assert payload['total_count']==record['total_count']==len(hits) and len({h['identifier'] for h in hits})==len(hits)
        entities.update(h['identifier'] for h in hits)
        dispositions.append(dict(sequence_id=sid,hits=len(hits),response_sha256=sha(responsepath),receipt_sha256=sha(rp)))
        print(number,'/ 39',sid,len(hits),'experimental entities',flush=True)
    assert len(dispositions)==39 and sha(fasta)==config['fasta_sha256'] and sha(sr)==config['source_receipt_sha256']
    (out/'sequence_dispositions.json').write_text(json.dumps(dispositions,indent=2)+'\n')
    (out/'polymer_entity_ids.txt').write_text('\n'.join(sorted(entities))+'\n')
    result=dict(status='complete_case_experimental_sequence_search_pending_readback',config_sha256=sha(cp),sequences=39,sequences_with_hits=sum(r['hits']>0 for r in dispositions),sequence_entity_hits=sum(r['hits'] for r in dispositions),unique_experimental_entities=len(entities),artifacts={name:sha(out/name) for name in ['sequence_dispositions.json','polymer_entity_ids.txt']},scope=config['scope'])
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
