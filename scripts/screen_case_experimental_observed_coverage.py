#!/usr/bin/env python3
"""Screen experimentally observed CA coverage for every case, chain and model."""
import csv,json,gzip,subprocess,time
import psutil
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
    launch_path=Path('metadata/case_experimental_ca_readback_launch_20260927.json')
    launch=json.loads(launch_path.read_text());sources[str(launch_path)]=sha(launch_path)
    while True:
        try:
            process=psutil.Process(launch['pid'])
            if process.create_time()!=launch['created'] or process.status()==psutil.STATUS_ZOMBIE:break
            assert process.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        print('Waiting for exact full CA readback',launch['pid'],flush=True);time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    assert sha(launch['cmdline'][1])==launch['script_sha256']
    mapping=base/'experimental_structures/whole-domain-case-ca-mapping-20260927-v1'
    audit_path=base/'experimental_structures/whole-domain-case-ca-readback-20260927-v1/receipt.json'
    ca_audit=json.loads(audit_path.read_text());assert ca_audit['status']=='complete_full_case_CA_raw_atom_and_grid_readback'
    sources[str(audit_path)]=sha(audit_path)
    mr=json.loads((mapping/'receipt.json').read_text());mc=json.loads((mapping/'config.json').read_text())
    assert mr['config_sha256']==sha(mapping/'config.json')
    assert ca_audit['source_hashes'][str(mapping/'receipt.json')]==sha(mapping/'receipt.json')
    sources[str(mapping/'receipt.json')]=sha(mapping/'receipt.json');sources[str(mapping/'config.json')]=sha(mapping/'config.json')
    observed=defaultdict(dict)
    for entry in mr['entry_receipt_sha256']:
        path=mapping/(entry+'.residues.tsv.gz')
        assert ca_audit['source_hashes'][str(path)]==sha(path);sources[str(path)]=sha(path)
        with gzip.open(path,'rt') as f:
            for row in csv.DictReader(f,delimiter='\t'):
                entity=entry+'_'+row['entity_id'];key=(row['model_number'],row['label_asym_id'])
                observed[entity].setdefault(key,set())
                if row['CA_status']=='unambiguous_full_occupancy_CA':observed[entity][key].add(int(row['label_seq_id']))
    out=base/'experimental_structures/whole-domain-case-observed-coverage-20260927-v1';out.mkdir(exist_ok=False)
    screens=[(n,c) for n in [30,50] for c in [50,70,90]];coverage={};hits=0;projected=0;excluded_hits=0;row_count=0
    output_path=out/'observed_alignment_interval_screens.tsv'
    with output_path.open('w') as handle:
        writer=None
        for d in disposition:
            p=search/(d['sequence_id']+'.json');assert sha(p)==d['response_sha256'];sources[str(p)]=sha(p)
            for hit in json.loads(p.read_text()).get('result_set',[]):
                entity=hit['identifier'];entry=entity.rsplit('_',1)[0]
                contexts=[c for service in hit['services'] for node in service['nodes'] for c in node['match_context']]
                for index,context in enumerate(contexts):
                    hits+=1
                    paired=paired_positions(context)
                    if entry in mc['deferred_entries']:
                        excluded_hits+=1;continue
                    assert entity in observed
                    for (model,chain),positions in observed[entity].items():
                        selected={q for q,s in paired if s in positions};projected+=1
                        for iid in seqintervals[d['sequence_id']]:
                            r=intervals[iid];domain=set(range(r['start'],r['end']+1))
                            inside=len(selected&domain);outside=len(selected-domain)
                            length=len(domain);remaining=r['original_length']-length
                            for ck,role in intervalcases[iid]:
                                for n,c in screens:
                                    anchor=inside>=n and inside*100>=c*length
                                    rest=outside>=n and remaining>0 and outside*100>=c*remaining
                                    screen=f'n{n}_c{c}'
                                    row=dict(zip(CASE,ck),role=role,interval_id=iid,sequence_id=d['sequence_id'],entity_id=entity,model_number=model,label_asym_id=chain,context_index=index,screen=screen,observed_domain_residues=inside,domain_length=length,observed_outside_residues=outside,outside_length=remaining,domain_pass=int(anchor),domain_and_outside_pass=int(anchor and rest))
                                    if writer is None:writer=csv.DictWriter(handle,list(row),delimiter='\t',lineterminator='\n');writer.writeheader()
                                    writer.writerow(row);row_count+=1
                                    key=ck,screen,entity,model,chain
                                    coverage.setdefault(key,{'domain':set(),'both':set()})
                                    if anchor:coverage[key]['domain'].add((role,iid))
                                    if anchor and rest:coverage[key]['both'].add((role,iid))
    assert hits==ar['alignment_contexts']==3091
    candidate_rows=[];summary=[]
    for r in dossiers:
        ck=tuple(r[k] for k in CASE);assert expected[ck]
        for n,c in screens:
            screen=f'n{n}_c{c}';domain_models=set();both_models=set()
            for (k,s,e,m,ch),values in sorted(coverage.items()):
                if k!=ck or s!=screen:continue
                domain=values['domain']==expected[ck];both=values['both']==expected[ck]
                if domain:domain_models.add((e,m,ch))
                if both:both_models.add((e,m,ch))
                candidate_rows.append(dict(zip(CASE,ck),screen=screen,entity_id=e,model_number=m,label_asym_id=ch,required_role_intervals=len(expected[ck]),passing_domain_role_intervals=len(values['domain']),passing_both_role_intervals=len(values['both']),shared_domain_pass=int(domain),shared_domain_and_outside_pass=int(both)))
            summary.append(dict(zip(CASE,ck),species_name=r['species_name'],screen=screen,shared_domain_chain_models=len(domain_models),shared_both_chain_models=len(both_models),shared_domain_entities=len({x[0] for x in domain_models}),shared_both_entities=len({x[0] for x in both_models}),shared_domain_entries=len({x[0].rsplit('_',1)[0] for x in domain_models}),shared_both_entries=len({x[0].rsplit('_',1)[0] for x in both_models})))
    assert len(summary)==78
    artifacts={output_path.name:sha(output_path)}
    for name,data in [('shared_chain_model_screens.tsv',candidate_rows),('case_summary.tsv',summary)]:
        p=out/name
        with p.open('w') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
        assert rows(p)==[{k:str(v) for k,v in r.items()} for r in data];artifacts[name]=sha(p)
    with output_path.open() as f:assert sum(1 for _ in csv.DictReader(f,delimiter='\t'))==row_count
    for path,digest in sources.items():assert sha(path)==digest
    result=dict(status='complete_case_experimental_observed_coverage',source_hashes=sources,script_sha256=sha(__file__),alignment_contexts=hits,excluded_integrative_alignment_contexts=excluded_hits,projected_chain_model_contexts=projected,interval_screen_rows=row_count,case_screen_rows=len(summary),artifacts=artifacts,scope='Only unambiguous standard-monomer full-occupancy CA positions projected through each query/subject alignment. A single deposited chain/model must pass every A/B/reference interval boundary; no merging chains, models or alignment contexts. Coverage eligibility only, not geometry, accuracy, function, training independence or independent experimental replication.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


def paired_positions(context):
    query,subject=context['query_beg']-1,context['subject_beg']-1
    pairs=[]
    qa,sa=context['query_aligned_seq'],context['subject_aligned_seq']
    assert len(qa)==len(sa)
    for q,s in zip(qa,sa):
        assert q!='-' or s!='-'
        if q!='-':query+=1
        if s!='-':subject+=1
        if q!='-' and s!='-':pairs.append((query,subject))
    assert query==context['query_end'] and subject==context['subject_end']
    assert len({q for q,s in pairs})==len(pairs)==len({s for q,s in pairs})
    return pairs


if __name__=='__main__':main()
