#!/usr/bin/env python3
"""Reconstruct original context gates and reference-pool direction states with SQL."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import itertools
from itertools import zip_longest
import json
from pathlib import Path
import sqlite3
from full_triad_context_contrast_sources import load_sources,MASKS,CORES,GATES,POLICIES,PHYSICAL_DIRECTIONS,STATES,SUMMARY_FIELDS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def run(path,output):
    plan=json.loads(Path(path).read_text());contexts,measurements,upstream,bindings=load_sources(plan,path)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text());rh=sha(rp)
    assert receipt['status']=='complete_full_triad_context_contrasts_pending_independent_readback' and receipt['plan_sha256']==sha(path)
    assert receipt['scientific_eligibility'] is False and receipt['scope']==plan['scope']
    assert receipt['source_gate_order']==GATES and receipt['policy_order']==POLICIES and receipt['physical_direction_order']==PHYSICAL_DIRECTIONS and receipt['direction_state_order']==STATES
    assert json.loads((root/'configuration.json').read_text())==dict(plan_sha256=sha(path),source_hashes=bindings)
    screens=[s['id'] for s in plan['screens']];db=sqlite3.connect(':memory:')
    db.execute('CREATE TABLE physical(t TEXT,m TEXT,c TEXT,models TEXT,'+','.join('q'+str(i)+' TEXT' for i in range(len(screens)))+',PRIMARY KEY(t,m,c))')
    db.execute('CREATE TABLE contexts(g TEXT,n INT,d TEXT,parent INT,lex_gene INT,lex_model INT,PRIMARY KEY(g,n,d))')
    db.execute('CREATE TABLE refs(d TEXT,m TEXT,c TEXT,s TEXT,lex INT,q0 TEXT,q1 TEXT,q2 TEXT)')
    groups=0;n=ties=0;guides=Counter();presence=Counter();counts=Counter()
    try:
        with gzip.open(measurements,'rt') as f:
            batch=[]
            for line in f:
                r=json.loads(line);m=r['mask_scenario'];c=r['core_scenario'];assert m in MASKS and c in CORES and r['role_order']==['a','b','reference']
                values=[]
                for sid in screens:
                    s=r['screens'][sid];assert type(s['all_selected_orders_pass']) is bool
                    if s['all_selected_orders_pass']:
                        assert s['qualified_direction']==r['numerical_tolerance_direction'] and s['qualified_direction'] in PHYSICAL_DIRECTIONS;values.append(s['qualified_direction'])
                    else:assert s['qualified_direction']=='excluded_by_quality';values.append(None)
                batch.append((r['triad_id'],m,c,json.dumps(r['models']),*values));groups+=1
                if len(batch)==5000:db.executemany('INSERT INTO physical VALUES('+','.join('?' for _ in batch[0])+')',batch);batch.clear()
            if batch:db.executemany('INSERT INTO physical VALUES('+','.join('?' for _ in batch[0])+')',batch)
        measured=db.execute('SELECT COUNT(DISTINCT t) FROM physical').fetchone()[0]
        assert measured==plan['expected']['measured_triads'] and groups==plan['expected']['measured_contrast_groups']
        assert db.execute('SELECT COUNT(*) FROM (SELECT t FROM physical GROUP BY t HAVING COUNT(*)!=9 OR COUNT(DISTINCT m)!=3 OR COUNT(DISTINCT c)!=3 OR COUNT(DISTINCT models)!=1)').fetchone()[0]==0
        number_chunks=(plan['expected']['target_contexts']+plan['checkpoint_contexts']-1)//plan['checkpoint_contexts']
        paths=[root/'checkpoints'/f'{i:06d}.jsonl.gz' for i in range(number_chunks)];assert set(paths)==set((root/'checkpoints').glob('*.jsonl.gz'))
        def checkpoints():
            for i,p in enumerate(paths):
                count=0
                with gzip.open(p,'rt') as f:
                    for line in f:count+=1;yield json.loads(line)
                assert count==min(plan['checkpoint_contexts'],plan['expected']['target_contexts']-i*plan['checkpoint_contexts'])
        canonical=lambda v:json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False)
        with gzip.open(contexts,'rt') as source,gzip.open(root/'context_contrasts.jsonl.gz','rt') as exported:
            for original_line,exported_line,checkpoint in zip_longest(source,exported,checkpoints()):
                assert original_line is not None and exported_line is not None and checkpoint is not None
                original=json.loads(original_line);actual=json.loads(exported_line)
                assert set(actual)=={'source_work_design','contrast_designs'} and canonical(actual['source_work_design'])==canonical(original)
                assert canonical(actual)==canonical(checkpoint) and set(actual['contrast_designs'])=={'availability','sequence_first'}
                work=original['source_design'];native=work['native_context'];guide=native['source_guide'];parent=native['parent_context_eligible'];assert type(parent) is bool
                db.execute('DELETE FROM refs');pending=[];first_gates={};reference_totals={}
                for design in ['availability','sequence_first']:
                    observed_design=actual['contrast_designs'][design]
                    assert set(observed_design)=={'references','context_policy_direction_states','context_policy_quality_flags','qualified_reference_counts','native_both_eligible_tie_direction_counts'}
                    refs=work['measurement_designs'][design];dispositions=original['triad_designs'][design];assert len(refs)==len(dispositions)==len(observed_design['references'])
                    lexical=refs[0]['reference'] if refs else None;reference_totals[design]=len(refs);first_gates[design]=None
                    db.execute('INSERT INTO contexts VALUES(?,?,?,?,?,?)',(guide,native['source_row_number'],design,int(parent),int(bool(lexical)),int(bool(lexical and lexical['reference_model']))))
                    for ix,(ref,disposition,observed) in enumerate(zip(refs,dispositions,observed_design['references'])):
                        chosen=ref['reference'];assert chosen['lexical_choice']==(ix==0) and disposition['reference_gene']==chosen['reference_gene']
                        models=[[work['duplicate_models'][role]['model_id'],work['duplicate_models'][role]['version']] for role in ['a','b']]
                        if chosen['reference_model']:
                            models.append([chosen['reference_model'],int(chosen['reference_version'])]);key=hashlib.sha256(json.dumps(models,separators=(',',':')).encode()).hexdigest()
                            distinct=len({tuple(v) for v in models})==3;ids=len({v[0] for v in models})==3
                            pairs=original['duplicate_pair_design']['measurement_disposition']=='in_full_primary_measurement_design' and all(ref['side_work'][s]['measurement_disposition']=='in_full_reference_measurement_design' for s in ['a','b'])
                        else:key=None;distinct=ids=pairs=False
                        assert disposition['triad_id']==key and disposition['three_distinct_versioned_models']==distinct and disposition['three_distinct_model_ids']==ids and disposition['all_three_pairs_in_measurement_design']==pairs
                        ready=bool(key and distinct and ids and pairs and parent and work['duplicate_comparison_status']=='queued_distinct_models')
                        own=bool(ready and chosen['native_coorthology'][guide]=='both');both=bool(ready and chosen['native_coorthology']['profile']==chosen['native_coorthology']['mafft']=='both');gates=[ready,own,both]
                        assert [disposition[g] for g in GATES]==gates
                        if ix==0:first_gates[design]=gates
                        physical={(m,c):(json.loads(payload),values) for m,c,payload,*values in db.execute('SELECT m,c,models,'+','.join('q'+str(i) for i in range(len(screens)))+' FROM physical WHERE t=?',(key,))}
                        present=bool(physical);assert not physical or len(physical)==9
                        if present:assert all(payload==models for payload,_ in physical.values())
                        if ready:assert present
                        qualified={}
                        for m in MASKS:
                            qualified[m]={}
                            for c in CORES:
                                qualified[m][c]={}
                                for i,sid in enumerate(screens):
                                    d=physical[m,c][1][i] if present else None;values=[d if gate else None for gate in gates]
                                    qualified[m][c][sid]=values;pending.append((design,m,c,sid,int(ix==0),*values))
                        expected=dict(reference_gene=chosen['reference_gene'],triad_id=key,physical_measurement_present=present,
                            physical_contrast_keys=[[key,m,c] for m,c in itertools.product(MASKS,CORES)] if present else None,
                            source_eligibility_flags=gates,qualified_direction_states=qualified)
                        assert canonical(observed)==canonical(expected)
                        presence[f'{guide}|{design}|parent={int(parent)}|measured={int(present)}|source_ready={int(ready)}']+=1;ties+=1
                db.executemany('INSERT INTO refs VALUES(?,?,?,?,?,?,?,?)',pending)
                sql='SELECT d,m,c,s,COUNT(*),COUNT(q0),COUNT(q1),COUNT(q2),MAX(CASE WHEN lex THEN q0 END),MAX(CASE WHEN lex THEN q1 END),MAX(CASE WHEN lex THEN q2 END),MIN(q2),MAX(q2),'+','.join('SUM(q2=?)' for _ in PHYSICAL_DIRECTIONS)+' FROM refs GROUP BY d,m,c,s'
                aggregates={(d,m,c,s):rest for d,m,c,s,*rest in db.execute(sql,PHYSICAL_DIRECTIONS)}
                for design in ['availability','sequence_first']:
                    states={};flags={};totals={};support={};num_refs=reference_totals[design]
                    for m in MASKS:
                        states[m]={};flags[m]={};totals[m]={};support[m]={}
                        for c in CORES:
                            states[m][c]={};flags[m][c]={};totals[m][c]={};support[m][c]={}
                            for sid in screens:
                                a=aggregates.get((design,m,c,sid),[0,0,0,0,None,None,None,None,None,0,0,0,0]);total,n0,n1,n2,l0,l1,l2,low,high,*directions=a
                                assert total==num_refs
                                pooled=low if low==high else 'reference_direction_disagreement'
                                lexical=[value if value is not None else 'no_tied_reference' if not num_refs else 'source_gate_excluded' if not first_gates[design][i] else 'geometry_quality_excluded' for i,value in enumerate([l0,l1,l2])]
                                any_state=pooled if n2 else 'no_eligible_reference' if num_refs else 'no_tied_reference'
                                all_state=pooled if num_refs and n2==num_refs else 'not_all_ties_eligible' if num_refs else 'no_tied_reference'
                                cell=lexical+[any_state,all_state];eligible=[v is not None for v in [l0,l1,l2]]+[n2>0,bool(num_refs and n2==num_refs)]
                                states[m][c][sid]=cell;flags[m][c][sid]=eligible;totals[m][c][sid]=[n0,n1,n2];support[m][c][sid]=[int(v or 0) for v in directions]
                                for policy,state in zip(POLICIES,cell):counts[guide,design,m,c,sid,policy,state]+=1
                    wanted=dict(context_policy_direction_states=states,context_policy_quality_flags=flags,qualified_reference_counts=totals,native_both_eligible_tie_direction_counts=support)
                    assert canonical({k:actual['contrast_designs'][design][k] for k in wanted})==canonical(wanted)
                n+=1;guides[guide]+=1
                if n%25000==0:print('Independent full context contrast SQL',n,'/',plan['expected']['target_contexts'],flush=True)
        assert n==plan['expected']['target_contexts'] and ties==plan['expected']['reference_tie_records'] and dict(guides)==upstream['guide_contexts']
        assert db.execute('SELECT COUNT(*) FROM contexts').fetchone()[0]==2*n
        for g,count in guides.items():
            for d in ['availability','sequence_first']:assert db.execute('SELECT MIN(n),MAX(n),COUNT(*) FROM contexts WHERE g=? AND d=?',(g,d)).fetchone()==(1,count,count)
        baseline={(g,d):[total,parent,gene,model] for g,d,total,parent,gene,model in db.execute('SELECT g,d,COUNT(*),SUM(parent),SUM(parent AND lex_gene),SUM(parent AND lex_model) FROM contexts GROUP BY g,d')}
        fields=['guide','design','mask_scenario','core_scenario','screen','policy','direction_state','context_count','source_contexts','parent_eligible_contexts','parent_eligible_lexical_gene_present','parent_eligible_lexical_model_present'];wanted_rows=[]
        for g,d,m,c,s,policy in itertools.product(sorted(guides),['availability','sequence_first'],MASKS,CORES,screens,POLICIES):
            assert sum(counts[g,d,m,c,s,policy,state] for state in STATES)==guides[g]
            for state in STATES:wanted_rows.append(dict(zip(fields,map(str,[g,d,m,c,s,policy,state,counts[g,d,m,c,s,policy,state],*baseline[g,d]]))))
        with (root/'context_contrast_counts.tsv').open() as f:
            reader=csv.DictReader(f,delimiter='\t');assert reader.fieldnames==fields and list(reader)==wanted_rows
    finally:db.close()
    cells=len(MASKS)*len(CORES)*len(screens)
    summary=dict(target_contexts=n,context_design_records=2*n,reference_tie_records=ties,duplicate_reference_links=2*ties,
        logical_reference_screen_decisions=ties*cells,context_screen_states=n*2*cells,context_policy_decisions=n*2*cells*5,
        summary_rows=len(wanted_rows),guide_contexts=dict(guides),reference_measurement_presence_counts=dict(presence),measured_triads=measured,measured_contrast_groups=groups,
        direction_counts={'|'.join(k):v for k,v in counts.items()})
    assert set(summary)==set(SUMMARY_FIELDS) and all(receipt[k]==v for k,v in summary.items())
    assert all(summary[k]==v for k,v in plan['expected'].items())
    bind(bindings,root/'configuration.json');assert receipt['source_hashes']==bindings
    expected_artifacts={'context_contrasts.jsonl.gz','context_contrast_counts.tsv',*[str(p.relative_to(root)) for p in paths]};assert set(receipt['artifacts'])==expected_artifacts
    for name,digest in receipt['artifacts'].items():bind(bindings,root/name,digest)
    bind(bindings,rp,rh);verify(bindings)
    result=dict(status='passed_full_triad_context_contrasts_original_gate_sql_readback',**summary,plan_sha256=sha(path),producer_receipt_sha256=rh,
        source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','direction_counts']}),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.plan,a.output)
