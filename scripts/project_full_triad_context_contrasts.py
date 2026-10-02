#!/usr/bin/env python3
"""Preserve original references and gates when assigning contrast direction states."""
import argparse
from collections import Counter
import csv
import fcntl
import gzip
import itertools
import json
from pathlib import Path
import shutil
from full_triad_context_contrast_sources import load_sources,MASKS,CORES,GATES,POLICIES,PHYSICAL_DIRECTIONS,STATES
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def physical_index(path,plan):
    index={};n=0
    with gzip.open(path,'rt') as f:
        for line in f:
            r=json.loads(line);key=(r['mask_scenario'],r['core_scenario']);assert key in set(itertools.product(MASKS,CORES))
            state=index.setdefault(r['triad_id'],dict(models=r['models'],leaves={}))
            assert state['models']==r['models'] and key not in state['leaves'] and r['role_order']==['a','b','reference']
            directions={}
            for s in plan['screens']:
                sid=s['id'];screen=r['screens'][sid];passed=screen['all_selected_orders_pass'];assert type(passed) is bool
                d=screen['qualified_direction'];assert d==(r['numerical_tolerance_direction'] if passed else 'excluded_by_quality')
                if passed:assert d in PHYSICAL_DIRECTIONS
                directions[sid]=d if passed else None
            state['leaves'][key]=directions;n+=1
    assert len(index)==plan['expected']['measured_triads'] and n==plan['expected']['measured_contrast_groups']
    assert all(set(v['leaves'])==set(itertools.product(MASKS,CORES)) for v in index.values())
    return index,n


def pool_direction(values):
    assert values and all(v in PHYSICAL_DIRECTIONS for v in values)
    return values[0] if len(set(values))==1 else 'reference_direction_disagreement'


def policy_cell(qualified,refs):
    totals=[sum(q[i] is not None for q in qualified) for i in range(3)]
    both=[q[2] for q in qualified if q[2] is not None];support=[both.count(d) for d in PHYSICAL_DIRECTIONS]
    lexical=[]
    for i in range(3):
        lexical.append('no_tied_reference' if not refs else qualified[0][i] if qualified[0][i] is not None else
                       'source_gate_excluded' if not refs[0]['source_eligibility_flags'][i] else 'geometry_quality_excluded')
    any_state=pool_direction(both) if both else 'no_eligible_reference' if refs else 'no_tied_reference'
    all_state=pool_direction(both) if refs and len(both)==len(refs) else 'not_all_ties_eligible' if refs else 'no_tied_reference'
    flags=[bool(qualified and qualified[0][i] is not None) for i in range(3)]+[bool(both),bool(refs and len(both)==len(refs))]
    return lexical+[any_state,all_state],flags,totals,support


def projection(original,index,plan,counts,presence,baseline):
    work=original['source_design'];native=work['native_context'];guide=native['source_guide'];parent=native['parent_context_eligible'];assert type(parent) is bool
    designs={};ties=0
    for design in ['availability','sequence_first']:
        refs=work['measurement_designs'][design];dispositions=original['triad_designs'][design];assert len(refs)==len(dispositions)
        lexical=refs[0]['reference'] if refs else None
        for label,value in [('source_contexts',1),('parent_eligible_contexts',parent),('parent_eligible_lexical_gene_present',bool(parent and lexical)),('parent_eligible_lexical_model_present',bool(parent and lexical and lexical['reference_model']))]:baseline[guide,design,label]+=value
        linked=[]
        for ix,(ref,disposition) in enumerate(zip(refs,dispositions)):
            chosen=ref['reference'];assert chosen['lexical_choice']==(ix==0) and chosen['reference_gene']==disposition['reference_gene']
            gates=[disposition[g] for g in GATES];assert all(type(g) is bool for g in gates) and (not any(gates) or parent)
            state=index.get(disposition['triad_id']);present=state is not None
            if present:
                models=[[work['duplicate_models'][role]['model_id'],work['duplicate_models'][role]['version']] for role in ['a','b']]+[[chosen['reference_model'],int(chosen['reference_version'])]]
                assert state['models']==models
            if gates[0]:assert present
            q={m:{c:{s['id']:[state['leaves'][m,c][s['id']] if present and gate else None for gate in gates] for s in plan['screens']} for c in CORES} for m in MASKS}
            linked.append(dict(reference_gene=chosen['reference_gene'],triad_id=disposition['triad_id'],physical_measurement_present=present,
                physical_contrast_keys=[[disposition['triad_id'],m,c] for m,c in itertools.product(MASKS,CORES)] if present else None,
                source_eligibility_flags=gates,qualified_direction_states=q))
            presence[f'{guide}|{design}|parent={int(parent)}|measured={int(present)}|source_ready={int(gates[0])}']+=1;ties+=1
        states={};flags={};qualified_counts={};support={}
        for m in MASKS:
            states[m]={};flags[m]={};qualified_counts[m]={};support[m]={}
            for c in CORES:
                states[m][c]={};flags[m][c]={};qualified_counts[m][c]={};support[m][c]={}
                for s in plan['screens']:
                    sid=s['id'];values=[r['qualified_direction_states'][m][c][sid] for r in linked]
                    cell,eligible,totals,directions=policy_cell(values,linked)
                    states[m][c][sid]=cell;flags[m][c][sid]=eligible;qualified_counts[m][c][sid]=totals;support[m][c][sid]=directions
                    for policy,d in zip(POLICIES,cell):counts[guide,design,m,c,sid,policy,d]+=1
        designs[design]=dict(references=linked,context_policy_direction_states=states,context_policy_quality_flags=flags,
            qualified_reference_counts=qualified_counts,native_both_eligible_tie_direction_counts=support)
    return dict(source_work_design=original,contrast_designs=designs),ties,guide


def execute(plan,path,contexts,measurements,upstream,bindings,out):
    index,groups=physical_index(measurements,plan);configuration=out/'configuration.json';config=dict(plan_sha256=sha(path),source_hashes=bindings)
    if configuration.exists():assert json.loads(configuration.read_text())==config
    else:
        temp=out/'configuration.tmp';temp.write_text(json.dumps(config,indent=2)+'\n');temp.replace(configuration)
    chunks=out/'checkpoints';chunks.mkdir(exist_ok=True);expected_chunks=set();artifacts={};chunk=[];number=0;reused=0
    def save():
        nonlocal number,reused
        if not chunk:return
        p=chunks/f'{number:06d}.jsonl.gz';expected_chunks.add(p.name);data=gzip.compress(('\n'.join(chunk)+'\n').encode(),compresslevel=1,mtime=0)
        if p.exists():assert p.read_bytes()==data;reused+=1
        else:
            temp=p.with_suffix('.tmp');temp.write_bytes(data);temp.replace(p)
        artifacts[str(p.relative_to(out))]=sha(p);chunk.clear();number+=1
    counts=Counter();presence=Counter();baseline=Counter();guides=Counter();n=ties=0
    with gzip.open(contexts,'rt') as source,gzip.open(out/'context_contrasts.jsonl.gz','wt',compresslevel=1) as dest:
        for line in source:
            r,nt,guide=projection(json.loads(line),index,plan,counts,presence,baseline);encoded=json.dumps(r,separators=(',',':'))
            dest.write(encoded+'\n');chunk.append(encoded);n+=1;ties+=nt;guides[guide]+=1
            if n%plan['checkpoint_contexts']==0:
                save();temp=out/'state.tmp';temp.write_text(json.dumps(dict(status='running_full_context_contrast_projection',completed_contexts=n,total_contexts=plan['expected']['target_contexts']))+'\n');temp.replace(out/'state.json')
                assert sum(p.stat().st_size for p in chunks.glob('*.jsonl.gz'))+dest.tell()<=plan['resources']['maximum_output_gib']*2**30
                assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
            if n%25000==0:print('Full original contrast context linkage',n,'/',plan['expected']['target_contexts'],flush=True)
        save()
    assert {p.name for p in chunks.glob('*.jsonl.gz')}==expected_chunks
    assert n==plan['expected']['target_contexts'] and ties==plan['expected']['reference_tie_records'] and dict(guides)==upstream['guide_contexts']
    fields=['guide','design','mask_scenario','core_scenario','screen','policy','direction_state','context_count','source_contexts','parent_eligible_contexts','parent_eligible_lexical_gene_present','parent_eligible_lexical_model_present'];rows=0
    with (out/'context_contrast_counts.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for g,d,m,c,s,policy in itertools.product(sorted(guides),['availability','sequence_first'],MASKS,CORES,[s['id'] for s in plan['screens']],POLICIES):
            assert sum(counts[g,d,m,c,s,policy,state] for state in STATES)==guides[g]
            for state in STATES:
                values=[g,d,m,c,s,policy,state,counts[g,d,m,c,s,policy,state],*[baseline[g,d,label] for label in fields[8:]]]
                w.writerow(dict(zip(fields,values)));rows+=1
    cells=len(MASKS)*len(CORES)*len(plan['screens'])
    summary=dict(target_contexts=n,context_design_records=2*n,reference_tie_records=ties,duplicate_reference_links=2*ties,
        logical_reference_screen_decisions=ties*cells,context_screen_states=n*2*cells,context_policy_decisions=n*2*cells*5,
        summary_rows=rows,guide_contexts=dict(guides),reference_measurement_presence_counts=dict(presence),measured_triads=len(index),measured_contrast_groups=groups,
        direction_counts={'|'.join(k):v for k,v in counts.items()})
    assert all(summary[k]==v for k,v in plan['expected'].items());bind(bindings,configuration);verify(bindings)
    for name in ['context_contrasts.jsonl.gz','context_contrast_counts.tsv']:artifacts[name]=sha(out/name)
    result=dict(status='complete_full_triad_context_contrasts_pending_independent_readback',**summary,plan_sha256=sha(path),checked_reused_chunks=reused,
        source_gate_order=GATES,policy_order=POLICIES,physical_direction_order=PHYSICAL_DIRECTIONS,direction_state_order=STATES,
        source_hashes=bindings,artifacts=artifacts,scientific_eligibility=False,scope=plan['scope'])
    assert not (out/'receipt.json').exists()
    temp=out/'receipt.tmp';temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(out/'receipt.json')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts','direction_counts']}),flush=True);return result


def run(path):
    plan=json.loads(Path(path).read_text());assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    contexts,measurements,upstream,bindings=load_sources(plan,path);out=Path(plan['output']);out.mkdir(exist_ok=True,parents=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (out/'receipt.json').exists()
        return execute(plan,path,contexts,measurements,upstream,bindings,out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
