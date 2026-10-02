#!/usr/bin/env python3
"""Original-context gate, joint-scenario and interrupted recovery contracts."""
import argparse
import contextlib
import copy
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys

import check_full_triad_context_geometry_cases as source_fixture
import project_full_triad_context_correspondence as producer
import readback_full_triad_context_correspondence as reader
from full_triad_context_correspondence_sources import FULL_BITS,MASKS,METHODS,DEFINITIONS,GATES
from run_ortholog_pair_guide_comparison import sha


def digest(models):return hashlib.sha256(json.dumps(models,separators=(',',':')).encode()).hexdigest()


def cases(original_plan):
    old=json.loads(Path(original_plan).read_text());temp=Path(original_plan).parent/'correspondence_contracts';temp.mkdir()
    context_path,old_physical,_,_=source_fixture.reader.load_sources(old,original_plan)
    with gzip.open(context_path,'rt') as f:original=[json.loads(s) for s in f]
    records=[]
    for guide in ['profile','mafft']:
        base=[r for r in original if r['source_design']['native_context']['source_guide']==guide]
        records.extend(base)
        for i in range(4):
            row=copy.deepcopy(base[0]);row['source_design']['native_context']['source_row_number']=11+i
            for design in ['availability','sequence_first']:
                chosen=row['source_design']['measurement_designs'][design][0]['reference'];chosen['reference_model']='r'+str(3+i)
                disposition=row['triad_designs'][design][0];disposition['triad_id']=digest([['a',1],['b',1],[chosen['reference_model'],1]])
            records.append(row)
    contexts=temp/'joint_contexts.jsonl.gz'
    with gzip.open(contexts,'wt') as f:
        for row in records:f.write(json.dumps(row)+'\n')
    physical=temp/'joint_comparisons.jsonl.gz';repeat=lambda bits:sum(bits<<(8*i) for i in range(6))
    with gzip.open(physical,'wt') as f:
        for i in range(7):
            models=[['a',1],['b',1],['r'+str(i),1]]
            for mask in MASKS[:2]:
                for method in METHODS:
                    for core in DEFINITIONS:
                        if i==0:bits=FULL_BITS
                        elif i==1:bits=FULL_BITS if mask=='full' else 0
                        elif i==2:bits=repeat(170 if mask=='full' else 85)
                        elif i==3:bits=FULL_BITS if method=='mafft_auto' else 0
                        elif i==4:bits=FULL_BITS if core=='reference_common' else 0
                        elif i==5:bits=repeat(170 if method=='famsa_default' else 85)
                        else:bits=FULL_BITS if (mask=='full' and method=='famsa_default') or (mask=='plddt70' and method=='mafft_auto') else 0
                        f.write(json.dumps(dict(triad_id=digest(models),models=models,mask=mask,sequence_method=method,mapping_definition=core,
                            screens={'n50_c70':dict(pair_pass_bits=bits,passing_pairs=bits.bit_count(),all_pairs_pass=bits==FULL_BITS,any_pair_pass=bits!=0)}))+'\n')
    ties=sum(len(r['triad_designs'][d]) for r in records for d in ['availability','sequence_first'])
    def sources(config,path):return contexts,physical,dict(guide_contexts=dict(profile=14,mafft=14)),{str(path):sha(path),str(contexts):sha(contexts),str(physical):sha(physical)}
    producer.load_sources=reader.load_sources=sources
    plan=dict(output=str(temp/'joint_projection'),screens=old['screens'],checkpoint_contexts=4,
        resources=dict(minimum_free_disk_gib=0,maximum_output_gib=1),pins={},
        expected=dict(target_contexts=28,reference_tie_records=ties,duplicate_reference_links=2*ties,measured_triads=7,measured_comparison_groups=56,
            logical_reference_screen_decisions=ties*27,context_screen_states=28*2*27,context_policy_decisions=28*2*27*5,summary_rows=540),
        scope='Synthetic full-context/source-gate software check, not a biological pilot. Original source proof I/O stubbed. Complete two-guide/two-design/27mask-method-core scenario grid, 48-bit leaves, missing/parent/native gates, lexical/any/all-tie policies, independent SHA/source reconstruction, complement-union bitmap and SQL summaries/checkpoint exports exercised.')
    pp=temp/'joint_plan.json';pp.write_text(json.dumps(plan))
    actual_projection=producer.projection;calls=0
    def interrupted(*args,**kw):
        nonlocal calls
        calls+=1
        if calls==6:raise RuntimeError('owned_context_fixture_interruption')
        return actual_projection(*args,**kw)
    producer.projection=interrupted
    try:producer.run(pp)
    except RuntimeError as e:assert str(e)=='owned_context_fixture_interruption'
    else:raise AssertionError('Fixture interruption absent')
    producer.projection=actual_projection
    made=producer.run(pp);checked=reader.run(pp,temp/'joint_readback.json')
    assert made['checked_reused_chunks']==1 and made['target_contexts']==checked['target_contexts']==28 and made['summary_rows']==540
    root=Path(plan['output'])
    with gzip.open(root/'context_correspondence.jsonl.gz','rt') as f:exported=[json.loads(s) for s in f]
    def flags(i,mask='both_masks',method='both_methods',core='both_cores'):
        return exported[i]['correspondence_designs']['availability']['context_policy_flags'][mask][method][core]['n50_c70']
    assert flags(0)==[True]*5 and flags(1)==[False]*5
    assert exported[1]['correspondence_designs']['availability']['references'][0]['physical_comparison_present'] is True
    assert flags(2)==[False,False,False,True,False] and flags(3)==[False]*5 and flags(4)==[False]*5
    assert flags(5)==[True,True,False,False,False] and flags(19)==[True,False,False,False,False]
    assert flags(6)==[True,True,True,True,False] and flags(6,'full')==[True]*5
    assert flags(7)==flags(8)==flags(9)==[False]*5 and flags(9,'full')==[True]*5
    assert flags(10,method='mafft_auto')==[True]*5 and flags(10,method='famsa_default')==flags(10)==[False]*5
    assert flags(11,core='reference_common')==[True]*5 and flags(11,core='cycle_consistent')==flags(11)==[False]*5
    ref=lambda g,i,ix=0:g[i]['correspondence_designs']['availability']['references'][ix]
    assert ref(exported,12)['pair_pass_bits']['full']['both_methods']['reference_common']['n50_c70']==0
    assert flags(13,'full','famsa_default')==[True]*5 and flags(13,'plddt70','mafft_auto')==[True]*5 and flags(13)==[False]*5
    receipt=json.loads((root/'receipt.json').read_text())
    mutations={
        'promoted_excluded_parent':lambda g,r:ref(g,1)['source_eligibility_flags'].__setitem__(0,True),
        'unmeasured_bitmap_as_zero':lambda g,r:ref(g,2)['pair_pass_bits']['full']['mafft_auto']['reference_common'].__setitem__('n50_c70',0),
        'favorable_reference_reselection':lambda g,r:g[2]['correspondence_designs']['availability']['references'].reverse(),
        'dropped_unmodeled_lexical_tie':lambda g,r:g[2]['correspondence_designs']['availability']['references'].pop(0),
        'vacuous_empty_all_ties':lambda g,r:g[3]['correspondence_designs']['availability']['context_policy_flags']['both_masks']['both_methods']['both_cores']['n50_c70'].__setitem__(4,True),
        'promoted_native_both_guides':lambda g,r:ref(g,5)['source_eligibility_flags'].__setitem__(2,True),
        'invented_both_mask_bitmap':lambda g,r:ref(g,8)['pair_pass_bits']['both_masks']['both_methods']['both_cores'].__setitem__('n50_c70',FULL_BITS),
        'invented_both_method_bitmap':lambda g,r:ref(g,10)['pair_pass_bits']['both_masks']['both_methods']['both_cores'].__setitem__('n50_c70',FULL_BITS),
        'invented_both_core_bitmap':lambda g,r:ref(g,11)['pair_pass_bits']['both_masks']['both_methods']['both_cores'].__setitem__('n50_c70',FULL_BITS),
        'union_complementary_method_orders':lambda g,r:ref(g,12)['pair_pass_bits']['full']['both_methods']['reference_common'].__setitem__('n50_c70',FULL_BITS),
        'favorable_mask_method_diagonal':lambda g,r:ref(g,13)['pair_pass_bits']['both_masks']['both_methods']['both_cores'].__setitem__('n50_c70',FULL_BITS),
        'changed_physical_key':lambda g,r:ref(g,0)['physical_comparison_keys'][0].__setitem__(0,'wrong'),
        'changed_source_taxon_context':lambda g,r:g[0]['source_work_design']['source_design']['native_context']['source'].update(taxon_id='changed'),
        'bitmap_boolean_type':lambda g,r:ref(g,11)['pair_pass_bits']['full']['mafft_auto']['cycle_consistent'].__setitem__('n50_c70',False),
        'qualified_integer_type':lambda g,r:ref(g,0)['qualified_all_pair_flags']['full']['mafft_auto']['reference_common']['n50_c70'].__setitem__(0,1),
        'duplicated_context':lambda g,r:g.insert(1,copy.deepcopy(g[0])),
        'missing_context':lambda g,r:g.pop(),
        'invented_scientific_eligibility':lambda g,r:r.update(scientific_eligibility=True),
    }
    rejected=[]
    for name,mutate in mutations.items():
        candidate=copy.deepcopy(exported);r=copy.deepcopy(receipt);mutate(candidate,r)
        case=root.parent/name;case.mkdir();cp=root.parent/(name+'_plan.json');config={**plan,'output':str(case)};cp.write_text(json.dumps(config))
        with gzip.open(case/'context_correspondence.jsonl.gz','wt') as f:
            for row in candidate:f.write(json.dumps(row)+'\n')
        (case/'context_correspondence_counts.tsv').write_bytes((root/'context_correspondence_counts.tsv').read_bytes())
        checkpoints=case/'checkpoints';checkpoints.mkdir()
        for i,start in enumerate(range(0,len(candidate),4)):
            (checkpoints/f'{i:06d}.jsonl.gz').write_bytes(gzip.compress(('\n'.join(json.dumps(row) for row in candidate[start:start+4])+'\n').encode(),mtime=0))
        source_hashes=sources(config,cp)[-1];configuration=case/'configuration.json';configuration.write_text(json.dumps(dict(plan_sha256=sha(cp),source_hashes=source_hashes)))
        r.update(plan_sha256=sha(cp),source_hashes={**source_hashes,str(configuration):sha(configuration)},artifacts={str(p.relative_to(case)):sha(p) for p in [case/'context_correspondence.jsonl.gz',case/'context_correspondence_counts.tsv',*checkpoints.glob('*.jsonl.gz')]})
        (case/'receipt.json').write_text(json.dumps(r))
        try:reader.run(cp,case/'readback.json')
        except AssertionError:rejected.append(name)
        else:raise AssertionError('Accepted false full-context export '+name)
    return dict(status='passed_full_triad_context_correspondence_software_contracts',target_contexts=28,reference_tie_records=ties,measured_triads=7,measured_comparison_groups=56,
        logical_reference_screen_decisions=made['logical_reference_screen_decisions'],context_screen_states=made['context_screen_states'],context_policy_decisions=made['context_policy_decisions'],summary_rows=540,
        checked_reused_chunks=1,rejected_rehashed_false_exports=rejected,scope=plan['scope'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    original=source_fixture.reader.run;result=None
    def hook(plan,output):
        nonlocal result
        checked=original(plan,output)
        if Path(output).name=='baseline-readback.json':result=cases(plan)
        return checked
    source_fixture.reader.run=hook;previous=sys.argv
    try:
        sys.argv=['check_full_triad_context_geometry_cases.py','--output',str(args.output.with_suffix('.source_fixture.json'))]
        with contextlib.redirect_stdout(io.StringIO()):source_fixture.main()
    finally:sys.argv=previous
    assert result is not None
    with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
