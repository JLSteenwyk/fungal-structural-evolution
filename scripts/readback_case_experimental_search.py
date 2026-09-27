#!/usr/bin/env python3
"""Verify all cached searches and reconstruct query alignment coverage for every hit context."""
import csv,hashlib,json,subprocess,time
from pathlib import Path
import psutil
from Bio import SeqIO
from screen_duplication_domain_alignment_coverage import sha


def main():
    lp=Path('metadata/case_experimental_sequence_search_launch_20260927.json');launch=json.loads(lp.read_text());lh=sha(lp)
    while True:
        try:
            p=psutil.Process(launch['pid'])
            if p.create_time()!=launch['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        print('Waiting for exact sequence search',launch['pid'],flush=True);time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    assert sha(lp)==lh and sha('scripts/search_case_experimental_homologs.py')==launch['script_sha256']
    root=Path('results/experimental_structures/whole-domain-case-sequence-search-20260927-v1');rp=root/'receipt.json';r=json.loads(rp.read_text());cp=root/'config.json';config=json.loads(cp.read_text());assert sha(cp)==r['config_sha256']
    fasta=Path('results/structural_comparisons/case-independent-control-inputs-20260927-v1/missing_independent_sequences.faa');assert sha(fasta)==config['fasta_sha256']
    sequences={x.id:str(x.seq) for x in SeqIO.parse(fasta,'fasta')};assert len(sequences)==39
    for name,digest in r['artifacts'].items():assert sha(root/name)==digest
    dispositions=json.loads((root/'sequence_dispositions.json').read_text());assert len(dispositions)==39 and {x['sequence_id'] for x in dispositions}==set(sequences)
    output=[];all_entities=set();hits_total=0;sources={str(x):sha(x) for x in [lp,rp,cp,fasta]}
    for d in dispositions:
        sid=d['sequence_id'];bp=root/(sid+'.json');brp=root/(sid+'.receipt.json');b=json.loads(brp.read_text());response=json.loads(bp.read_text())
        assert sha(brp)==d['receipt_sha256'] and sha(bp)==d['response_sha256']==b['response_sha256'] and b['config_sha256']==sha(cp)
        sources[str(bp)]=sha(bp);sources[str(brp)]=sha(brp)
        query=b['query'];parameters=query['query']['parameters']
        assert parameters==dict(sequence_type='protein',value=sequences[sid],identity_cutoff=.3,evalue_cutoff=1e-5)
        assert query['request_options']['results_content_type']==['experimental'] and query['request_options']['return_all_hits'] is True
        hits=response.get('result_set',[]);assert len(hits)==response['total_count']==d['hits']==b['total_count']
        assert len({x['identifier'] for x in hits})==len(hits);hits_total+=len(hits)
        for hit in hits:
            entity=hit['identifier'];all_entities.add(entity);contexts=[]
            for service in hit['services']:
                assert service['service_type']=='sequence'
                for node in service['nodes']:contexts.extend(node['match_context'])
            assert contexts
            for index,context in enumerate(contexts):
                qa=context['query_aligned_seq'];sa=context['subject_aligned_seq'];qb=context['query_beg'];qe=context['query_end'];sb=context['subject_beg'];se=context['subject_end']
                assert len(qa)==len(sa)==context['alignment_length'] and not any(a==b=='-' for a,b in zip(qa,sa))
                assert 1<=qb<=qe<=len(sequences[sid])==context['query_length']
                assert qa.replace('-','')==sequences[sid][qb-1:qe]
                assert len(sa.replace('-',''))==se-sb+1 and 1<=sb<=se<=context['subject_length']
                paired=sum(a!='-' and b!='-' for a,b in zip(qa,sa));identical=sum(a==b and a!='-' for a,b in zip(qa,sa))
                output.append(dict(sequence_id=sid,entity_id=entity,context_index=index,query_begin=qb,query_end=qe,subject_begin=sb,subject_end=se,query_length=len(sequences[sid]),subject_length=context['subject_length'],aligned_columns=len(qa),paired_residues=paired,identical_residues=identical,query_span_coverage=(qe-qb+1)/len(sequences[sid]),query_paired_coverage=paired/len(sequences[sid]),subject_span_coverage=(se-sb+1)/context['subject_length'],identity_per_aligned_column=identical/len(qa),identity_per_paired_residue=identical/paired,reported_identity=context['sequence_identity'],reported_evalue=context['evalue']))
    assert hits_total==r['sequence_entity_hits'] and len(all_entities)==r['unique_experimental_entities'] and all_entities==set((root/'polymer_entity_ids.txt').read_text().split())
    out=Path('results/experimental_structures/whole-domain-case-sequence-readback-20260927-v1');out.mkdir(exist_ok=False);p=out/'alignment_contexts.tsv'
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]),delimiter='\t');w.writeheader();w.writerows(output)
    with p.open() as f:assert list(csv.DictReader(f,delimiter='\t'))==[{k:str(v) for k,v in row.items()} for row in output]
    for path,digest in sources.items():assert sha(path)==digest
    result=dict(status='passed_full_case_experimental_sequence_search_readback',terminal_state=state,source_hashes=sources,checker_sha256=sha(__file__),sequences=39,sequences_with_hits=r['sequences_with_hits'],sequence_entity_hits=hits_total,unique_entities=len(all_entities),alignment_contexts=len(output),artifacts={p.name:sha(p)},scope='All sequence queries, response counts and alignment contexts checked; query aligned residues match original sequence substrings exactly. Subject span lengths checked, but subject sequence identity awaits independent entity metadata. Coverage is alignment coverage, not experimentally observed coordinate coverage. No homology or functional claim is established solely by search hits.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
