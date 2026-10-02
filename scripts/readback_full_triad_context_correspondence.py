#!/usr/bin/env python3
"""Independently reconstruct original gates, scenario intersections and context policies."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import sqlite3

from full_triad_context_correspondence_sources import load_sources,MASKS,METHODS,METHOD_SCENARIOS,DEFINITIONS,CORE_SCENARIOS,FULL_BITS,GATES,POLICIES,SUMMARY_FIELDS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path,output):
    plan=json.loads(Path(plan_path).read_text());contexts,measurements,upstream,bindings=load_sources(plan,plan_path)
    root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    assert r['status']=='complete_full_triad_context_correspondence_pending_independent_readback' and r['plan_sha256']==sha(plan_path)
    assert r['masks']==MASKS and r['methods']==METHOD_SCENARIOS and r['cores']==CORE_SCENARIOS and r['source_gate_order']==GATES and r['context_policy_order']==POLICIES
    assert r['scope']==plan['scope'] and r['scientific_eligibility'] is False
    assert json.loads((root/'configuration.json').read_text())==dict(plan_sha256=sha(plan_path),source_hashes=bindings)
    screens=[s['id'] for s in plan['screens']];db=sqlite3.connect(':memory:')
    db.execute('CREATE TABLE physical(t TEXT,m TEXT,method TEXT,core TEXT,models TEXT,'+','.join('b'+str(i)+' INT' for i in range(len(screens)))+',PRIMARY KEY(t,m,method,core))')
    db.execute('CREATE TABLE contexts(g TEXT,n INT,d TEXT,parent INT,lex_gene INT,lex_model INT,PRIMARY KEY(g,n,d))')
    db.execute('CREATE TABLE refs(d TEXT,m TEXT,method TEXT,core TEXT,s TEXT,lex INT,ready INT,own INT,both INT,bits INT)')
    groups=0;expected_artifacts={'context_correspondence.jsonl.gz','context_correspondence_counts.tsv'}
    try:
        with gzip.open(measurements,'rt') as f:
            batch=[]
            for line in f:
                row=json.loads(line);assert row['mask'] in MASKS[:2] and row['sequence_method'] in METHODS and row['mapping_definition'] in DEFINITIONS
                bits=[row['screens'][sid]['pair_pass_bits'] for sid in screens];assert all(type(b) is int and 0<=b<=FULL_BITS for b in bits)
                for sid,b in zip(screens,bits):
                    screen=row['screens'][sid];assert screen['passing_pairs']==bin(b).count('1') and screen['all_pairs_pass'] is (b==FULL_BITS) and screen['any_pair_pass'] is (b!=0)
                batch.append([row['triad_id'],row['mask'],row['sequence_method'],row['mapping_definition'],json.dumps(row['models']),*bits]);groups+=1
                if len(batch)==10000:
                    db.executemany('INSERT INTO physical VALUES('+','.join('?' for _ in batch[0])+')',batch);batch.clear()
            if batch:db.executemany('INSERT INTO physical VALUES('+','.join('?' for _ in batch[0])+')',batch)
        measured=db.execute('SELECT COUNT(DISTINCT t) FROM physical').fetchone()[0]
        assert measured==plan['expected']['measured_triads'] and groups==plan['expected']['measured_comparison_groups']
        assert db.execute('SELECT COUNT(*) FROM (SELECT t FROM physical GROUP BY t HAVING COUNT(*)!=8 OR COUNT(DISTINCT m)!=2 OR COUNT(DISTINCT method)!=2 OR COUNT(DISTINCT core)!=2 OR COUNT(DISTINCT models)!=1)').fetchone()[0]==0
        expected_chunks=(plan['expected']['target_contexts']+plan['checkpoint_contexts']-1)//plan['checkpoint_contexts']
        chunk_paths=[root/'checkpoints'/f'{i:06d}.jsonl.gz' for i in range(expected_chunks)]
        assert set(root.joinpath('checkpoints').glob('*.jsonl.gz'))==set(chunk_paths)
        def archived_records():
            for i,path in enumerate(chunk_paths):
                expected_artifacts.add(str(path.relative_to(root)));n=0
                with gzip.open(path,'rt') as f:
                    for line in f:n+=1;yield json.loads(line)
                assert n==min(plan['checkpoint_contexts'],plan['expected']['target_contexts']-i*plan['checkpoint_contexts'])
        counts=Counter();guides=Counter();presence=Counter();n=ties=0
        with gzip.open(contexts,'rt') as source,gzip.open(root/'context_correspondence.jsonl.gz','rt') as exported:
            for left,right,checkpoint in zip_longest(source,exported,archived_records()):
                assert left is not None and right is not None and checkpoint is not None,'Missing or extra context/checkpoint'
                original=json.loads(left);actual=json.loads(right);assert actual==checkpoint
                assert set(actual)=={'source_work_design','correspondence_designs'} and actual['source_work_design']==original
                assert json.dumps(actual['source_work_design'],sort_keys=True,separators=(',',':'))==json.dumps(original,sort_keys=True,separators=(',',':'))
                assert set(actual['correspondence_designs'])=={'availability','sequence_first'}
                work=original['source_design'];native=work['native_context'];guide=native['source_guide'];parent=native['parent_context_eligible'];assert type(parent) is bool
                db.execute('DELETE FROM refs');pending=[]
                for design in ['availability','sequence_first']:
                    observed_design=actual['correspondence_designs'][design];assert set(observed_design)=={'references','context_policy_flags'}
                    refs=work['measurement_designs'][design];dispositions=original['triad_designs'][design]
                    assert len(refs)==len(dispositions)==len(observed_design['references'])
                    lexical=refs[0]['reference'] if refs else None
                    db.execute('INSERT INTO contexts VALUES(?,?,?,?,?,?)',(guide,native['source_row_number'],design,int(parent),int(bool(lexical)),int(bool(lexical and lexical['reference_model']))))
                    for ix,(ref,disposition,observed) in enumerate(zip(refs,dispositions,observed_design['references'])):
                        chosen=ref['reference'];assert chosen['lexical_choice']==(ix==0) and disposition['reference_gene']==chosen['reference_gene']
                        models=[[work['duplicate_models'][role]['model_id'],work['duplicate_models'][role]['version']] for role in ['a','b']]
                        if chosen['reference_model']:
                            models.append([chosen['reference_model'],int(chosen['reference_version'])])
                            key=hashlib.sha256(json.dumps(models,separators=(',',':')).encode()).hexdigest()
                            distinct_versioned=len({tuple(m) for m in models})==3;distinct_ids=len({m[0] for m in models})==3
                            all_pairs=original['duplicate_pair_design']['measurement_disposition']=='in_full_primary_measurement_design' and all(ref['side_work'][side]['measurement_disposition']=='in_full_reference_measurement_design' for side in ['a','b'])
                        else:key=None;distinct_versioned=distinct_ids=all_pairs=False
                        assert disposition['triad_id']==key and disposition['three_distinct_versioned_models']==distinct_versioned and disposition['three_distinct_model_ids']==distinct_ids and disposition['all_three_pairs_in_measurement_design']==all_pairs
                        ready=bool(key and distinct_versioned and distinct_ids and all_pairs and parent and work['duplicate_comparison_status']=='queued_distinct_models')
                        own=bool(ready and chosen['native_coorthology'][guide]=='both');both=bool(ready and chosen['native_coorthology']['profile']==chosen['native_coorthology']['mafft']=='both')
                        gates=[ready,own,both];assert [disposition[gate] for gate in GATES]==gates
                        physical={(m,method,core):(json.loads(payload),bits) for m,method,core,payload,*bits in db.execute('SELECT m,method,core,models,'+','.join('b'+str(i) for i in range(len(screens)))+' FROM physical WHERE t=?',(key,))}
                        present=bool(physical);assert not physical or len(physical)==8
                        if present:assert all(payload==models for payload,_ in physical.values())
                        if ready:assert present
                        bitmaps={};qualified={}
                        for mask in MASKS:
                            bitmaps[mask]={};qualified[mask]={}
                            for method in METHOD_SCENARIOS:
                                bitmaps[mask][method]={};qualified[mask][method]={}
                                for core in CORE_SCENARIOS:
                                    cells={};flags={}
                                    for i,sid in enumerate(screens):
                                        bits=None
                                        if present:
                                            values=[v[i] for (m,q,d),(_,v) in physical.items() if (mask=='both_masks' or mask==m) and (method=='both_methods' or method==q) and (core=='both_cores' or core==d)]
                                            assert len(values)==(2 if mask=='both_masks' else 1)*(2 if method=='both_methods' else 1)*(2 if core=='both_cores' else 1)
                                            missing=0
                                            for value in values:missing |=FULL_BITS^value
                                            bits=FULL_BITS^missing
                                        cells[sid]=bits;flags[sid]=[bool(bits==FULL_BITS and gate) for gate in gates]
                                        pending.append((design,mask,method,core,sid,int(ix==0),*map(int,gates),bits))
                                    bitmaps[mask][method][core]=cells;qualified[mask][method][core]=flags
                        expected=dict(reference_gene=chosen['reference_gene'],triad_id=key,physical_comparison_present=present,
                            physical_comparison_keys=[[key,m,q,d] for m in MASKS[:2] for q in METHODS for d in DEFINITIONS] if present else None,
                            source_eligibility_flags=gates,pair_pass_bits=bitmaps,qualified_all_pair_flags=qualified)
                        assert observed==expected and type(observed['physical_comparison_present']) is bool
                        assert all(v is None or type(v) is int and 0<=v<=FULL_BITS for m in MASKS for q in METHOD_SCENARIOS for d in CORE_SCENARIOS for v in observed['pair_pass_bits'][m][q][d].values())
                        assert all(type(v) is bool for v in observed['source_eligibility_flags'])
                        assert all(type(v) is bool for m in MASKS for q in METHOD_SCENARIOS for d in CORE_SCENARIOS for sid in screens for v in observed['qualified_all_pair_flags'][m][q][d][sid])
                        presence[f'{guide}|{design}|parent={int(parent)}|measured={int(present)}|source_ready={int(ready)}']+=1;ties+=1
                db.executemany('INSERT INTO refs VALUES(?,?,?,?,?,?,?,?,?,?)',pending)
                aggregates={}
                sql=f'SELECT d,m,method,core,s,COUNT(*),SUM(lex AND ready AND COALESCE(bits,0)={FULL_BITS}),SUM(lex AND own AND COALESCE(bits,0)={FULL_BITS}),SUM(lex AND both AND COALESCE(bits,0)={FULL_BITS}),SUM(both AND COALESCE(bits,0)={FULL_BITS}) FROM refs GROUP BY d,m,method,core,s'
                for design,mask,method,core,sid,count,lexical,own,both,passed in db.execute(sql):aggregates[design,mask,method,core,sid]=[bool(lexical),bool(own),bool(both),passed>0,count>0 and passed==count]
                for design in ['availability','sequence_first']:
                    policies={}
                    for mask in MASKS:
                        policies[mask]={}
                        for method in METHOD_SCENARIOS:
                            policies[mask][method]={}
                            for core in CORE_SCENARIOS:
                                cells={}
                                for sid in screens:
                                    flags=aggregates.get((design,mask,method,core,sid),[False]*5);cells[sid]=flags
                                    for policy,flag in zip(POLICIES,flags):counts[guide,design,mask,method,core,sid,policy]+=flag
                                policies[mask][method][core]=cells
                    observed=actual['correspondence_designs'][design]['context_policy_flags'];assert observed==policies
                    assert all(type(v) is bool for m in MASKS for q in METHOD_SCENARIOS for d in CORE_SCENARIOS for sid in screens for v in observed[m][q][d][sid])
                n+=1;guides[guide]+=1
                if n%25000==0:print('Independent SQL full context correspondence',n,'/',plan['expected']['target_contexts'],flush=True)
        assert n==plan['expected']['target_contexts'] and ties==plan['expected']['reference_tie_records'] and dict(guides)==upstream['guide_contexts']
        assert db.execute('SELECT COUNT(*) FROM contexts').fetchone()[0]==2*n
        for guide,count in guides.items():
            for design in ['availability','sequence_first']:assert db.execute('SELECT MIN(n),MAX(n),COUNT(*) FROM contexts WHERE g=? AND d=?',(guide,design)).fetchone()==(1,count,count)
        baseline={(g,d):(count,parent,genes,models) for g,d,count,parent,genes,models in db.execute('SELECT g,d,COUNT(*),SUM(parent),SUM(parent AND lex_gene),SUM(parent AND lex_model) FROM contexts GROUP BY g,d')}
        fields=['guide','design','mask','sequence_method_scenario','core_scenario','screen','policy','source_contexts','parent_eligible_contexts','parent_eligible_lexical_gene_present','parent_eligible_lexical_model_present','passed_contexts']
        expected_rows=[dict(zip(fields,map(str,[*key,*baseline[key[:2]],count]))) for key,count in sorted(counts.items())]
        with (root/'context_correspondence_counts.tsv').open() as f:reader=csv.DictReader(f,delimiter='\t');assert reader.fieldnames==fields and list(reader)==expected_rows
    finally:db.close()
    cells=len(MASKS)*len(METHOD_SCENARIOS)*len(CORE_SCENARIOS)*len(screens)
    summary=dict(target_contexts=n,context_design_records=2*n,reference_tie_records=ties,duplicate_reference_links=2*ties,
        logical_reference_screen_decisions=ties*cells,context_screen_states=n*2*cells,context_policy_decisions=n*2*cells*len(POLICIES),summary_rows=len(counts),
        guide_contexts=dict(guides),reference_measurement_presence_counts=dict(presence),measured_triads=measured,measured_comparison_groups=groups)
    assert set(summary)==set(SUMMARY_FIELDS) and all(r[k]==v for k,v in summary.items())
    assert all(summary[k]==v for k,v in plan['expected'].items())
    bind(bindings,root/'configuration.json');assert r['source_hashes']==bindings
    assert set(r['artifacts'])==expected_artifacts
    bind(bindings,rp,rh)
    for name,digest in r['artifacts'].items():bind(bindings,root/name,digest)
    verify(bindings)
    result=dict(status='passed_full_triad_context_correspondence_sql_readback',**summary,plan_sha256=sha(plan_path),producer_receipt_sha256=rh,
        source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.plan,a.output)
