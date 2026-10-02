#!/usr/bin/env python3
"""Original-context/reference-pool contracts and interrupted checkpoint recovery."""
import argparse
import contextlib
import copy
import gzip
import hashlib
import io
import itertools
import json
from pathlib import Path
import sys
import check_full_triad_context_geometry_cases as source_fixture
import project_full_triad_context_contrasts as producer
import readback_full_triad_context_contrasts as reader
from full_triad_context_contrast_sources import MASKS,CORES,GATES
from run_ortholog_pair_guide_comparison import sha


def digest(models):return hashlib.sha256(json.dumps(models,separators=(',',':')).encode()).hexdigest()


def contracts(old_path):
    old=json.loads(Path(old_path).read_text());temp=Path(old_path).parent/'contrast_contracts';temp.mkdir()
    contexts,_,_,_=source_fixture.reader.load_sources(old,old_path)
    with gzip.open(contexts,'rt') as f:base=[json.loads(line) for line in f]
    records=[]
    for guide in ['profile','mafft']:
        native=[r for r in base if r['source_design']['native_context']['source_guide']==guide];records.extend(native)
        for i in range(2):
            row=copy.deepcopy(native[0]);row['source_design']['native_context']['source_row_number']=11+i
            for design in ['availability','sequence_first']:
                chosen=row['source_design']['measurement_designs'][design][0]['reference'];chosen['reference_model']='r'+str(3+i)
                row['triad_designs'][design][0]['triad_id']=digest([['a',1],['b',1],[chosen['reference_model'],1]])
            records.append(row)
    source=temp/'contexts.jsonl.gz';physical=temp/'physical.jsonl.gz'
    with gzip.open(source,'wt') as f:
        for row in records:f.write(json.dumps(row)+'\n')
    with gzip.open(physical,'wt') as f:
        for i in range(5):
            models=[['a',1],['b',1],['r'+str(i),1]]
            for m,c in itertools.product(MASKS,CORES):
                d='positive' if i in [0,3] else 'negative' if i==1 else 'within_numerical_tolerance' if i==4 else 'positive' if m=='full' else 'negative' if m=='plddt70' and c=='reference_common' else 'positive' if m=='plddt70' and c=='cycle_consistent' else 'sign_uncertain'
                passed=i!=3
                f.write(json.dumps(dict(triad_id=digest(models),models=models,role_order=['a','b','reference'],mask_scenario=m,core_scenario=c,
                    numerical_tolerance_direction=d,screens={'n50_c70':dict(all_selected_orders_pass=passed,qualified_direction=d if passed else 'excluded_by_quality')}))+'\n')
    ties=sum(len(r['triad_designs'][d]) for r in records for d in ['availability','sequence_first']);assert ties==52
    def sources(config,path):return source,physical,dict(guide_contexts=dict(profile=12,mafft=12)),{str(path):sha(path),str(source):sha(source),str(physical):sha(physical)}
    producer.load_sources=reader.load_sources=sources
    plan=dict(output=str(temp/'baseline'),screens=old['screens'],checkpoint_contexts=4,resources=dict(minimum_free_disk_gib=0,maximum_output_gib=1),
        expected=dict(target_contexts=24,reference_tie_records=52,duplicate_reference_links=104,measured_triads=5,measured_contrast_groups=45,
            logical_reference_screen_decisions=468,context_screen_states=432,context_policy_decisions=2160,summary_rows=1800),
        scope='Synthetic full original context/reference-pool software contracts; production source-closure I/O stubbed, not biological evidence.')
    pp=temp/'plan.json';pp.write_text(json.dumps(plan));original=producer.projection;calls=0
    def interrupted(*a,**kw):
        nonlocal calls
        calls+=1
        if calls==6:raise RuntimeError('owned_context_contrast_interruption')
        return original(*a,**kw)
    producer.projection=interrupted
    try:producer.run(pp)
    except RuntimeError as e:assert str(e)=='owned_context_contrast_interruption'
    else:raise AssertionError('Missing owned interruption')
    finally:producer.projection=original
    receipt=producer.run(pp);checked=reader.run(pp,temp/'readback.json');assert receipt['checked_reused_chunks']==1 and checked['context_policy_decisions']==2160
    root=Path(plan['output'])
    with gzip.open(root/'context_contrasts.jsonl.gz','rt') as f:exported=[json.loads(line) for line in f]
    def design(g,i):return g[i]['contrast_designs']['availability']
    def states(g,i,m='both_masks',c='both_cores'):return design(g,i)['context_policy_direction_states'][m][c]['n50_c70']
    assert states(exported,0)==['positive']*5 and states(exported,1)==['source_gate_excluded']*3+['no_eligible_reference','not_all_ties_eligible']
    assert design(exported,1)['references'][0]['physical_measurement_present'] is True
    assert states(exported,2)==['source_gate_excluded']*3+['positive','not_all_ties_eligible'] and states(exported,3)==['no_tied_reference']*5
    assert states(exported,5)==['positive','positive','source_gate_excluded','no_eligible_reference','not_all_ties_eligible']
    assert states(exported,17)==['positive','source_gate_excluded','source_gate_excluded','no_eligible_reference','not_all_ties_eligible']
    assert states(exported,6)==['positive']*3+['reference_direction_disagreement']*2
    assert design(exported,6)['context_policy_quality_flags']['both_masks']['both_cores']['n50_c70']==[True]*5
    assert design(exported,6)['native_both_eligible_tie_direction_counts']['both_masks']['both_cores']['n50_c70']==[1,1,0,0]
    assert states(exported,8)==['sign_uncertain']*5 and states(exported,8,'full')==['positive']*5 and states(exported,9)==['negative']*5
    assert states(exported,10)==['geometry_quality_excluded']*3+['no_eligible_reference','not_all_ties_eligible'] and states(exported,11)==['within_numerical_tolerance']*5
    ref=lambda g,i:design(g,i)['references'][0]
    mutations={
        'promoted_excluded_parent':lambda g,r:ref(g,1)['source_eligibility_flags'].__setitem__(0,True),
        'missing_direction_as_positive':lambda g,r:ref(g,2)['qualified_direction_states']['both_masks']['both_cores']['n50_c70'].__setitem__(0,'positive'),
        'favorable_lexical_reselection':lambda g,r:design(g,2)['references'].reverse(),
        'dropped_missing_lexical_tie':lambda g,r:design(g,2)['references'].pop(0),
        'vacuous_empty_all_ties':lambda g,r:design(g,3)['context_policy_quality_flags']['both_masks']['both_cores']['n50_c70'].__setitem__(4,True),
        'promoted_native_both_gate':lambda g,r:ref(g,5)['source_eligibility_flags'].__setitem__(2,True),
        'reference_disagreement_as_positive':lambda g,r:states(g,6).__setitem__(4,'positive'),
        'partial_reference_pool_as_complete':lambda g,r:states(g,2).__setitem__(4,'positive'),
        'favorable_mask_direction':lambda g,r:states(g,8).__setitem__(2,'positive'),
        'nearzero_as_positive':lambda g,r:states(g,11).__setitem__(2,'positive'),
        'invented_quality_qualification':lambda g,r:design(g,10)['context_policy_quality_flags']['both_masks']['both_cores']['n50_c70'].__setitem__(2,True),
        'wrong_tie_direction_support':lambda g,r:design(g,6)['native_both_eligible_tie_direction_counts']['both_masks']['both_cores']['n50_c70'].__setitem__(1,0),
        'wrong_qualified_reference_count':lambda g,r:design(g,6)['qualified_reference_counts']['both_masks']['both_cores']['n50_c70'].__setitem__(2,1),
        'changed_physical_key':lambda g,r:ref(g,0)['physical_contrast_keys'][0].__setitem__(0,'wrong'),
        'changed_source_taxon':lambda g,r:g[0]['source_work_design']['source_design']['native_context']['source'].update(taxon_id='changed'),
        'integer_quality_flag':lambda g,r:design(g,0)['context_policy_quality_flags']['both_masks']['both_cores']['n50_c70'].__setitem__(0,1),
        'boolean_reference_count':lambda g,r:design(g,0)['qualified_reference_counts']['both_masks']['both_cores']['n50_c70'].__setitem__(0,True),
        'duplicate_context':lambda g,r:g.insert(0,copy.deepcopy(g[0])),
        'missing_context':lambda g,r:g.pop(),
        'invented_scientific_eligibility':lambda g,r:r.update(scientific_eligibility=True),
    }
    rejected=[]
    for label,mutate in mutations.items():
        case=temp/label;case.mkdir();config={**plan,'output':str(case)};cp=temp/(label+'.json');cp.write_text(json.dumps(config))
        candidate=copy.deepcopy(exported);r=copy.deepcopy(receipt);mutate(candidate,r)
        with gzip.open(case/'context_contrasts.jsonl.gz','wt') as f:
            for record in candidate:f.write(json.dumps(record)+'\n')
        chunks=case/'checkpoints';chunks.mkdir()
        for i,start in enumerate(range(0,len(candidate),4)):(chunks/f'{i:06d}.jsonl.gz').write_bytes(gzip.compress(('\n'.join(json.dumps(v) for v in candidate[start:start+4])+'\n').encode(),mtime=0))
        (case/'context_contrast_counts.tsv').write_bytes((root/'context_contrast_counts.tsv').read_bytes())
        bindings=sources(config,cp)[-1];configuration=case/'configuration.json';configuration.write_text(json.dumps(dict(plan_sha256=sha(cp),source_hashes=bindings)))
        r.update(plan_sha256=sha(cp),source_hashes={**bindings,str(configuration):sha(configuration)},artifacts={str(p.relative_to(case)):sha(p) for p in [case/'context_contrasts.jsonl.gz',case/'context_contrast_counts.tsv',*chunks.glob('*.jsonl.gz')]})
        (case/'receipt.json').write_text(json.dumps(r))
        try:reader.run(cp,case/'readback.json')
        except AssertionError:rejected.append(label)
        else:raise AssertionError('Accepted false full context contrast export '+label)
    return dict(status='passed_full_triad_context_contrast_software_contracts',**plan['expected'],checked_reused_chunks=1,rejected_rehashed_false_exports=rejected,scientific_eligibility=False,scope=plan['scope'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();original=source_fixture.reader.run;result=None
    def hook(plan,output):
        nonlocal result
        checked=original(plan,output)
        if Path(output).name=='baseline-readback.json':result=contracts(plan)
        return checked
    source_fixture.reader.run=hook;previous=sys.argv
    try:
        sys.argv=['check_full_triad_context_geometry_cases.py','--output',str(a.output.with_suffix('.source_fixture.json'))]
        with contextlib.redirect_stdout(io.StringIO()):source_fixture.main()
    finally:sys.argv=previous
    assert result is not None;a.output.parent.mkdir(exist_ok=True,parents=True)
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
