#!/usr/bin/env python3
"""Preserve every original context and tie under joint correspondence qualification."""
import argparse
from collections import Counter
import csv
import fcntl
import gzip
import itertools
import json
from pathlib import Path
import shutil

from full_triad_context_correspondence_sources import load_sources,leaf_selectors,MASKS,METHODS,METHOD_SCENARIOS,DEFINITIONS,CORE_SCENARIOS,FULL_BITS,GATES,POLICIES
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def physical_index(measurements,plan):
    index={};groups=0
    with gzip.open(measurements,'rt') as f:
        for line in f:
            row=json.loads(line);tid=row['triad_id'];key=(row['mask'],row['sequence_method'],row['mapping_definition'])
            assert key[0] in MASKS[:2] and key[1] in METHODS and key[2] in DEFINITIONS
            state=index.setdefault(tid,dict(models=row['models'],leaves={}))
            assert state['models']==row['models'] and key not in state['leaves']
            leaf={s['id']:row['screens'][s['id']]['pair_pass_bits'] for s in plan['screens']}
            assert all(type(v) is int and 0<=v<=FULL_BITS for v in leaf.values())
            for sid,bits in leaf.items():
                screen=row['screens'][sid]
                assert screen['passing_pairs']==bits.bit_count() and screen['all_pairs_pass'] is (bits==FULL_BITS) and screen['any_pair_pass'] is (bits!=0)
            state['leaves'][key]=leaf;groups+=1
    expected=set(itertools.product(MASKS[:2],METHODS,DEFINITIONS))
    assert len(index)==plan['expected']['measured_triads'] and groups==plan['expected']['measured_comparison_groups']
    assert all(set(s['leaves'])==expected for s in index.values())
    return index,groups


def projection(original,index,plan,presence,baseline,counts):
    native=original['source_design']['native_context'];guide=native['source_guide'];parent=native['parent_context_eligible']
    assert type(parent) is bool
    projected={};ties=0
    for design in ['availability','sequence_first']:
        references=original['triad_designs'][design];source_refs=original['source_design']['measurement_designs'][design]
        assert len(references)==len(source_refs);linked=[];lexical=source_refs[0]['reference'] if source_refs else None
        for label,value in [('source_contexts',1),('parent_eligible_contexts',parent),('parent_eligible_lexical_gene_present',bool(parent and lexical)),('parent_eligible_lexical_model_present',bool(parent and lexical and lexical['reference_model']))]:baseline[guide,design,label]+=value
        for ix,(ref,source_ref) in enumerate(zip(references,source_refs)):
            chosen=source_ref['reference'];assert chosen['lexical_choice']==(ix==0) and chosen['reference_gene']==ref['reference_gene']
            flags=[ref[gate] for gate in GATES];assert all(type(v) is bool for v in flags) and (not any(flags) or parent)
            state=index.get(ref['triad_id']);measured=state is not None
            if measured:
                models=[[original['source_design']['duplicate_models'][role]['model_id'],original['source_design']['duplicate_models'][role]['version']] for role in ['a','b']]
                models.append([chosen['reference_model'],int(chosen['reference_version'])]);assert state['models']==models
            if flags[0]:assert measured
            bits={};qualified={}
            for mask in MASKS:
                bits[mask]={};qualified[mask]={}
                for method in METHOD_SCENARIOS:
                    bits[mask][method]={};qualified[mask][method]={}
                    for core in CORE_SCENARIOS:
                        cells={};q={}
                        for screen in plan['screens']:
                            sid=screen['id'];value=None
                            if measured:
                                value=FULL_BITS
                                for key in itertools.product(*leaf_selectors(mask,method,core)):value &=state['leaves'][key][sid]
                            cells[sid]=value;q[sid]=[bool(flag and value==FULL_BITS) for flag in flags]
                        bits[mask][method][core]=cells;qualified[mask][method][core]=q
            keys=[[ref['triad_id'],mask,method,core] for mask in MASKS[:2] for method in METHODS for core in DEFINITIONS] if measured else None
            linked.append(dict(reference_gene=ref['reference_gene'],triad_id=ref['triad_id'],physical_comparison_present=measured,physical_comparison_keys=keys,
                source_eligibility_flags=flags,pair_pass_bits=bits,qualified_all_pair_flags=qualified))
            presence[f'{guide}|{design}|parent={int(parent)}|measured={int(measured)}|source_ready={int(flags[0])}']+=1;ties+=1
        policies={}
        for mask in MASKS:
            policies[mask]={}
            for method in METHOD_SCENARIOS:
                policies[mask][method]={}
                for core in CORE_SCENARIOS:
                    cells={}
                    for screen in plan['screens']:
                        sid=screen['id'];q=[r['qualified_all_pair_flags'][mask][method][core][sid] for r in linked];first=q[0] if q else [False]*3
                        cells[sid]=[*first,any(v[2] for v in q),bool(q and all(v[2] for v in q))]
                        for policy,value in zip(POLICIES,cells[sid]):counts[guide,design,mask,method,core,sid,policy]+=value
                    policies[mask][method][core]=cells
        projected[design]=dict(references=linked,context_policy_flags=policies)
    return dict(source_work_design=original,correspondence_designs=projected),ties,guide


def execute(plan,plan_path,contexts,measurements,upstream,bindings,out):
    index,groups=physical_index(measurements,plan)
    configuration=out/'configuration.json';value=dict(plan_sha256=sha(plan_path),source_hashes=bindings)
    if configuration.exists():assert json.loads(configuration.read_text())==value,'Changed context checkpoint sources'
    else:
        with configuration.open('x') as f:f.write(json.dumps(value,indent=2)+'\n')
    checkpoints=out/'checkpoints';checkpoints.mkdir(exist_ok=True);artifacts={};expected_files=set();chunk=[];chunk_number=0;reused=0
    checkpoint_bytes=sum(p.stat().st_size for p in checkpoints.glob('*.jsonl.gz'))
    def commit_chunk():
        nonlocal chunk_number,reused,checkpoint_bytes
        if not chunk:return
        name=f'{chunk_number:06d}.jsonl.gz';expected_files.add(name);path=checkpoints/name
        payload=gzip.compress(('\n'.join(chunk)+'\n').encode(),compresslevel=1,mtime=0)
        if path.exists():assert path.read_bytes()==payload,'False or changed context checkpoint';reused+=1
        else:
            temp=path.with_suffix('.tmp');temp.write_bytes(payload);temp.replace(path);checkpoint_bytes+=len(payload)
        artifacts[str(path.relative_to(out))]=sha(path);chunk.clear();chunk_number+=1
    guides=Counter();presence=Counter();baseline=Counter();counts=Counter();n=ties=0
    names=['context_correspondence.jsonl.gz','context_correspondence_counts.tsv']
    with gzip.open(contexts,'rt') as source,gzip.open(out/names[0],'wt',compresslevel=1) as dest:
        for line in source:
            record,nt,guide=projection(json.loads(line),index,plan,presence,baseline,counts)
            encoded=json.dumps(record,separators=(',',':'));dest.write(encoded+'\n');chunk.append(encoded);n+=1;ties+=nt;guides[guide]+=1
            if n%plan['checkpoint_contexts']==0:
                commit_chunk();state=dict(status='running_full_context_correspondence_projection',completed_contexts=n,total_contexts=plan['expected']['target_contexts'],checked_reused_chunks=reused)
                temp=out/'state.tmp';temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(out/'state.json')
                assert checkpoint_bytes+dest.tell()<=plan['resources']['maximum_output_gib']*2**30
                assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
            if n%25000==0:print('Full original context correspondence projection',n,'/',plan['expected']['target_contexts'],flush=True)
        commit_chunk()
    assert {p.name for p in checkpoints.glob('*.jsonl.gz')}==expected_files
    assert n==plan['expected']['target_contexts'] and ties==plan['expected']['reference_tie_records'] and dict(guides)==upstream['guide_contexts']
    fields=['guide','design','mask','sequence_method_scenario','core_scenario','screen','policy','source_contexts','parent_eligible_contexts','parent_eligible_lexical_gene_present','parent_eligible_lexical_model_present','passed_contexts']
    with (out/names[1]).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for key,count in sorted(counts.items()):
            value=dict(zip(fields[:7],key));value.update({label:baseline[key[0],key[1],label] for label in fields[7:-1]});value['passed_contexts']=count;w.writerow(value)
    cells=len(MASKS)*len(METHOD_SCENARIOS)*len(CORE_SCENARIOS)*len(plan['screens'])
    summary=dict(target_contexts=n,context_design_records=2*n,reference_tie_records=ties,duplicate_reference_links=2*ties,
        logical_reference_screen_decisions=ties*cells,context_screen_states=n*2*cells,context_policy_decisions=n*2*cells*len(POLICIES),summary_rows=len(counts),
        guide_contexts=dict(guides),reference_measurement_presence_counts=dict(presence),measured_triads=len(index),measured_comparison_groups=groups)
    assert all(summary[k]==v for k,v in plan['expected'].items())
    bind(bindings,configuration);verify(bindings)
    for name in names:artifacts[name]=sha(out/name)
    result=dict(status='complete_full_triad_context_correspondence_pending_independent_readback',**summary,plan_sha256=sha(plan_path),checked_reused_chunks=reused,
        masks=MASKS,methods=METHOD_SCENARIOS,cores=CORE_SCENARIOS,source_gate_order=GATES,context_policy_order=POLICIES,
        source_hashes=bindings,artifacts=artifacts,scientific_eligibility=False,scope=plan['scope'])
    with (out/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}),flush=True);return result


def run(plan_path):
    plan=json.loads(Path(plan_path).read_text());assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    contexts,measurements,upstream,bindings=load_sources(plan,plan_path);out=Path(plan['output']);out.mkdir(exist_ok=True,parents=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (out/'receipt.json').exists(),'Completed context output immutable'
        return execute(plan,plan_path,contexts,measurements,upstream,bindings,out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
