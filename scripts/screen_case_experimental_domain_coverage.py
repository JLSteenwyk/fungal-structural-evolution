#!/usr/bin/env python3
"""Screen paired query residues within/outside every case domain, retaining shared entities."""
import csv,json
from collections import defaultdict
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha

CASE=['family','gene_a','gene_b','pfam_accession']


def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    base=Path('results');sources={}
    def checked(root,name):
        rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/name;assert sha(p)==r['artifacts'][name]
        sources[str(rp)]=sha(rp);sources[str(p)]=sha(p);return p
    structural=base/'structural_comparisons';case=structural/'whole-domain-case-dossiers-20260927-v1'
    dossiers=rows(checked(case,'case_dossiers.tsv'));links=rows(checked(case,'domain_reference_links.tsv'))
    membership=defaultdict(set)
    for r in links:membership[r['triad_key']].add(tuple(r[k] for k in CASE))
    intervalcases=defaultdict(set)
    for line in checked(structural/'duplication-domain-common-residues-20260927-v1','triads.jsonl').open():
        r=json.loads(line)
        for ck in membership.get(r['triad_key'],[]):
            for role in ['a','b','reference']:intervalcases[r['interval_'+role]].add((ck,role))
    intervals={}
    for line in checked(structural/'duplication-domain-inputs-20260926-v1','inputs.jsonl').open():
        r=json.loads(line)
        if r['mask']=='full' and r['interval_id'] in intervalcases:intervals[r['interval_id']]=r
    modelseq={(m['model_id'],m['version']):'S'+m['sequence_sha256'] for m in json.loads(checked(structural/'case-independent-control-inputs-20260927-v1','model_provenance.json').read_text())}
    seqintervals=defaultdict(list);expected=defaultdict(set)
    for iid,r in intervals.items():
        sid=modelseq[r['model_id'],r['version']];seqintervals[sid].append(iid)
        for ck,role in intervalcases[iid]:expected[ck].add((role,iid))
    search=base/'experimental_structures/whole-domain-case-sequence-search-20260927-v1'
    disposition=json.loads(checked(search,'sequence_dispositions.json').read_text())
    audit=base/'experimental_structures/whole-domain-case-sequence-readback-20260927-v1/receipt.json'
    ar=json.loads(audit.read_text());assert ar['status']=='passed_full_case_experimental_sequence_search_readback';sources[str(audit)]=sha(audit)
    output=[];coverage={};screens=[(n,c) for n in [30,50] for c in [50,70,90]];hits=0
    for d in disposition:
        p=search/(d['sequence_id']+'.json');assert sha(p)==d['response_sha256'];sources[str(p)]=sha(p)
        for hit in json.loads(p.read_text()).get('result_set',[]):
            contexts=[c for s in hit['services'] for node in s['nodes'] for c in node['match_context']]
            for index,context in enumerate(contexts):
                hits+=1;pos=context['query_beg']-1;paired=[]
                for a,b in zip(context['query_aligned_seq'],context['subject_aligned_seq']):
                    if a!='-':pos+=1
                    if a!='-' and b!='-':paired.append(pos)
                assert pos==context['query_end'] and len(paired)==len(set(paired))
                for iid in seqintervals[d['sequence_id']]:
                    r=intervals[iid];inside=sum(r['start']<=p<=r['end'] for p in paired);outside=len(paired)-inside
                    length=r['end']-r['start']+1;remaining=r['original_length']-length
                    assert inside==len(set(paired)&set(range(r['start'],r['end']+1)))
                    for ck,role in intervalcases[iid]:
                        for n,c in screens:
                            anchor=inside>=n and inside*100>=c*length
                            rest=outside>=n and remaining>0 and outside*100>=c*remaining
                            screen=f'n{n}_c{c}'
                            output.append(dict(zip(CASE,ck),role=role,interval_id=iid,sequence_id=d['sequence_id'],entity_id=hit['identifier'],context_index=index,screen=screen,paired_domain_residues=inside,domain_length=length,paired_outside_residues=outside,outside_length=remaining,domain_pass=int(anchor),domain_and_outside_pass=int(anchor and rest)))
                            key=ck,screen,hit['identifier'];coverage.setdefault(key,{'domain':set(),'both':set()})
                            if anchor:coverage[key]['domain'].add((role,iid))
                            if anchor and rest:coverage[key]['both'].add((role,iid))
    assert hits==ar['alignment_contexts']==3091
    entities=[];summary=[]
    for r in dossiers:
        ck=tuple(r[k] for k in CASE);assert expected[ck]
        for n,c in screens:
            screen=f'n{n}_c{c}';counts={'domain':0,'both':0}
            for (k,s,e),values in sorted(coverage.items()):
                if k!=ck or s!=screen:continue
                domain=values['domain']==expected[ck];both=values['both']==expected[ck]
                counts['domain']+=domain;counts['both']+=both
                entities.append(dict(zip(CASE,ck),screen=screen,entity_id=e,required_role_intervals=len(expected[ck]),passing_domain_role_intervals=len(values['domain']),passing_both_role_intervals=len(values['both']),shared_domain_pass=int(domain),shared_domain_and_outside_pass=int(both)))
            summary.append(dict(zip(CASE,ck),species_name=r['species_name'],screen=screen,shared_domain_entities=counts['domain'],shared_domain_and_outside_entities=counts['both']))
    assert len(summary)==78
    out=base/'experimental_structures/whole-domain-case-domain-coverage-20260927-v1';out.mkdir(exist_ok=False);artifacts={}
    for name,data in [('alignment_interval_screens.tsv',output),('shared_entity_screens.tsv',entities),('case_summary.tsv',summary)]:
        p=out/name
        with p.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
        assert rows(p)==[{k:str(v) for k,v in r.items()} for r in data];artifacts[name]=sha(p)
    for path,digest in sources.items():assert sha(path)==digest
    result=dict(status='complete_case_experimental_domain_alignment_coverage',source_hashes=sources,script_sha256=sha(__file__),alignment_contexts=hits,alignment_interval_screen_rows=len(output),shared_entity_screen_rows=len(entities),case_screen_rows=len(summary),artifacts=artifacts,scope='Coverage counts only query residues paired to a subject residue; gaps do not contribute. A shared entity must pass every A/B/reference interval boundary for a case. All search hits retained. Multiple contexts, if present, qualify independently and are never unioned to manufacture coverage. Experimental observed-coordinate coverage and construct review remain outstanding.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
