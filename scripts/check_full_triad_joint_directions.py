#!/usr/bin/env python3
"""Full joint-grid raw-fit/bitmap/direction contracts and checkpoint recovery."""
import argparse
import copy
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import tempfile
import summarize_full_triad_joint_directions as producer
import readback_full_triad_joint_directions as reader
from compare_full_triad_correspondence import sequence_summary,pair,joint_summary
from summarize_full_triad_order_robustness import summarize
from full_triad_joint_direction_sources import MASKS,METHODS,DEFINITIONS,PERMUTATIONS,FULL_BITS
from full_triad_correspondence_comparison_sources import SEQUENCE_NUMERIC_FIELDS
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    with tempfile.TemporaryDirectory(prefix='fungal-full-joint-directions-') as folder:
        root=Path(folder);paths={name:root/(name+('.tsv.gz' if name.startswith('raw_') else '.jsonl.gz')) for name in ['raw_sequence','raw_structural','sequence_groups','structural_groups','comparison_groups']}
        screens=[dict(id='n50_c70'),dict(id='n30_c50')];config=dict(screens=screens);triads=[];raw=[[],[]];groups=[[],[],[]]
        for i in range(14):
            models=[[f'a{i}',1],[f'b{i}',2],[f'r{i}',3]]
            if i==1:models[0],models[1]=models[1],models[0]
            tid=hashlib.sha256(json.dumps(models,separators=(',',':')).encode()).hexdigest();triad=dict(triad_id=tid,models=models);triads.append(triad);seq={};struct={}
            for kind,axes in [(0,METHODS),(1,DEFINITIONS)]:
                for m,x in itertools.product(MASKS[:2],axes):
                    for order in range(6 if kind==0 else 8):
                        v=1.;status='computed_unique_at_numeric_tolerance'
                        if i==1:v=-1.
                        elif i==2 and kind==1:v=-1.
                        elif i==3:v=1. if m=='full' else -1.
                        elif i==4 and kind==0:v=1. if x=='famsa_default' else -1.
                        elif i==5 and kind==1:v=1. if x=='reference_common' else -1.
                        elif i==6 and kind==0 and order==5:v=-1.
                        elif i==7 and kind==1 and order==7:v=-1.
                        elif i==8 and kind==0 and m=='plddt70' and x=='mafft_auto' and order==5:v=None;status='source_excluded'
                        elif i==9 and kind==1:status='computed_nonunique_at_numeric_tolerance'
                        elif i==10:v=5e-10
                        elif i==11:v=1e-9
                        row=dict(triad_id=tid,mask=m,fit_status=status,triples_sha256='fixture-triples',source_exclusions='')
                        row.update({field:1. for field in SEQUENCE_NUMERIC_FIELDS});row.update(rmsd_ar=2.+v if v is not None else '',rmsd_br=2. if v is not None else '',rmsd_ar_minus_br=v if v is not None else '')
                        for role,model in zip(['a','b','reference'],models):row['model_'+role],row['version_'+role]=model
                        if kind==0:row.update(method=x,permutation=PERMUTATIONS[order],sequence_native_status='valid_sequence_alignment',native_payload_sha256='fixture-native')
                        else:
                            row['mapping_definition']=x
                            for edge,digit in zip(['ab','ar','br'],f'{order:03b}'):row['order_'+edge]=int(digit)
                        for s in screens:
                            sid=s['id'];core=int(status=='computed_unique_at_numeric_tolerance');inherited=1
                            if sid=='n50_c70' and kind==0:
                                if i==12:inherited=int(order%2==(0 if x=='famsa_default' else 1))
                                if i==13:inherited=int((m=='full' and x=='famsa_default') or (m=='plddt70' and x=='mafft_auto'))
                            row.update({sid+'_pass':core*inherited,sid+'_core_pass':core,sid+'_three_pair_pass':inherited,sid+'_exclusions':'' if core*inherited else 'fixture_excluded'})
                        raw[kind].append(row);value={k:str(v) for k,v in row.items()}
                        if kind==0:seq[m,x,PERMUTATIONS[order]]=value
                        else:struct[m,x,order]=value
            for m,q in itertools.product(MASKS[:2],METHODS):groups[0].append(sequence_summary(triad,m,q,[seq[m,q,p] for p in PERMUTATIONS],config))
            for m,c in itertools.product(MASKS[:2],DEFINITIONS):groups[1].append(dict(triad_id=tid,models=models,mask=m,mapping_definition=c,role_order=['a','b','reference'],order_bits_layout=[list(o) for o in itertools.product([0,1],repeat=3)],**summarize([struct[m,c,i] for i in range(8)],screens)))
            for m,q,c in itertools.product(MASKS[:2],METHODS,DEFINITIONS):
                pairs=[pair(triad,m,q,p,c,i,seq[m,q,p],struct[m,c,i],[(0,0,0)],[(0,0,0)],config) for p in PERMUTATIONS for i in range(8)]
                groups[2].append(joint_summary(triad,m,q,c,pairs,config))
        for kind,name in [(0,'raw_sequence'),(1,'raw_structural')]:
            with gzip.open(paths[name],'wt') as f:
                w=csv.DictWriter(f,fieldnames=list(raw[kind][0]),delimiter='\t');w.writeheader();w.writerows(reversed(raw[kind]))
        for name,records in zip(['sequence_groups','structural_groups','comparison_groups'],groups):
            with gzip.open(paths[name],'wt') as f:
                for record in records:f.write(json.dumps(record)+'\n')
        plan=dict(output=str(root/'baseline'),screens=screens,checkpoint_triads=2,contrast_numerical_tolerance_angstrom=1e-9,resources=dict(minimum_free_disk_gib=0,maximum_output_gib=1),
            expected=dict(measured_triads=14,sequence_fit_rows=336,structural_fit_rows=448,source_sequence_groups=56,source_structural_groups=56,source_comparison_groups=112,joint_direction_groups=378,screen_decisions=756,summary_rows=378),
            scope='Synthetic full27-scenario joint direction software grid; production source closure I/O stubbed, not biological evidence.')
        pp=root/'plan.json';pp.write_text(json.dumps(plan))
        def sources(config,path):return paths,triads,{str(p):sha(p) for p in [path,*paths.values()]}
        producer.load_sources=reader.load_sources=sources;original=producer.scenario;calls=0
        def interrupted(*a,**kw):
            nonlocal calls
            calls+=1
            if calls==56:raise RuntimeError('owned_joint_direction_interruption')
            return original(*a,**kw)
        producer.scenario=interrupted
        try:producer.run(pp)
        except RuntimeError as e:assert str(e)=='owned_joint_direction_interruption'
        else:raise AssertionError('Missing owned interruption')
        finally:producer.scenario=original
        receipt=producer.run(pp);checked=reader.run(pp,root/'readback.json');assert receipt['checked_reused_chunks']==1 and checked['screen_decisions']==756
        try:producer.run(pp)
        except AssertionError:pass
        else:raise AssertionError('Accepted completed producer restart')
        out=Path(plan['output'])
        with gzip.open(out/'joint_directions.jsonl.gz','rt') as f:exported=[json.loads(line) for line in f]
        joint=lambda i:exported[i*27+26]
        assert [joint(i)['joint_envelope']['numerical_tolerance_direction'] for i in range(14)]==['positive','negative','sign_uncertain','sign_uncertain','sign_uncertain','sign_uncertain','sign_uncertain','sign_uncertain','unavailable','nonunique_fit','within_numerical_tolerance','within_numerical_tolerance','positive','positive']
        assert joint(0)['joint_envelope']['expected_orders']==56 and joint(8)['joint_envelope']['available_orders']==55
        assert joint(2)['source_direction_category_relation']=='different_categories' and joint(3)['source_direction_category_relation']=='matching_categories'
        assert joint(11)['joint_envelope']['strict_zero_direction']=='positive'
        assert joint(12)['screens']['n50_c70']['joint_pair_pass_bits']==joint(13)['screens']['n50_c70']['joint_pair_pass_bits']==0
        assert joint(12)['screens']['n30_c50']['joint_pair_pass_bits']==FULL_BITS
        g=lambda records,i:records[i*27+26]
        mutations={
            'favorable_source_selection':lambda v,r:g(v,2)['joint_envelope'].update(numerical_tolerance_direction='positive'),
            'matching_uncertainty_as_stable':lambda v,r:g(v,3)['joint_envelope'].update(numerical_tolerance_direction='positive'),
            'favorable_method':lambda v,r:g(v,4)['sequence_envelope'].update(numerical_tolerance_direction='positive'),
            'favorable_core':lambda v,r:g(v,5)['structural_envelope'].update(numerical_tolerance_direction='positive'),
            'dropped_sequence_order':lambda v,r:g(v,6)['sequence_envelope'].update(minimum_contrast_angstrom=1.),
            'dropped_structural_order':lambda v,r:g(v,7)['structural_envelope'].update(minimum_contrast_angstrom=1.),
            'incomplete_as_complete':lambda v,r:g(v,8)['joint_envelope'].update(available_orders=56),
            'nonunique_as_unique':lambda v,r:g(v,9)['joint_envelope'].update(all_selected_orders_unique=True),
            'nearzero_as_positive':lambda v,r:g(v,10)['joint_envelope'].update(numerical_tolerance_direction='positive'),
            'union_complementary_methods':lambda v,r:g(v,12)['screens']['n50_c70'].update(joint_pair_pass_bits=FULL_BITS),
            'favorable_mask_method_diagonal':lambda v,r:g(v,13)['screens']['n50_c70'].update(joint_pair_pass_bits=FULL_BITS),
            'invented_quality':lambda v,r:g(v,12)['screens']['n50_c70'].update(all_selected_pairs_pass=True,qualified_joint_direction='positive'),
            'changed_models':lambda v,r:g(v,0)['models'].reverse(),
            'wrong_source_key':lambda v,r:g(v,0)['sequence_group_keys'][0].__setitem__(0,'wrong'),
            'boolean_count':lambda v,r:g(v,0)['joint_envelope'].update(available_orders=True),
            'integer_unique':lambda v,r:g(v,0)['joint_envelope'].update(all_selected_orders_unique=1),
            'missing_group':lambda v,r:v.pop(),
            'duplicate_group':lambda v,r:v.insert(0,copy.deepcopy(v[0])),
            'invented_scientific_eligibility':lambda v,r:r.update(scientific_eligibility=True),
            'changed_tolerance':lambda v,r:r.update(contrast_numerical_tolerance_angstrom=0.),
        }
        rejected=[]
        for name,mutate in mutations.items():
            case=root/name;case.mkdir();config={**plan,'output':str(case)};cp=root/(name+'.json');cp.write_text(json.dumps(config));candidate=copy.deepcopy(exported);r=copy.deepcopy(receipt);mutate(candidate,r)
            with gzip.open(case/'joint_directions.jsonl.gz','wt') as f:
                for record in candidate:f.write(json.dumps(record)+'\n')
            chunks=case/'checkpoints';chunks.mkdir()
            for i,start in enumerate(range(0,len(candidate),54)):(chunks/f'{i:06d}.jsonl.gz').write_bytes(gzip.compress(('\n'.join(json.dumps(v) for v in candidate[start:start+54])+'\n').encode(),mtime=0))
            (case/'joint_direction_counts.tsv').write_bytes((out/'joint_direction_counts.tsv').read_bytes());bindings=sources(config,cp)[-1];configuration=case/'configuration.json';configuration.write_text(json.dumps(dict(plan_sha256=sha(cp),source_hashes=bindings)))
            r.update(plan_sha256=sha(cp),source_hashes={**bindings,str(configuration):sha(configuration)},artifacts={str(p.relative_to(case)):sha(p) for p in [case/'joint_directions.jsonl.gz',case/'joint_direction_counts.tsv',*chunks.glob('*.jsonl.gz')]})
            (case/'receipt.json').write_text(json.dumps(r))
            try:reader.run(cp,case/'readback.json')
            except AssertionError:rejected.append(name)
            else:raise AssertionError('Accepted false joint direction export '+name)
        result=dict(status='passed_full_triad_joint_direction_software_contracts',**plan['expected'],checked_reused_chunks=1,completed_restart_refused=True,
            raw_fit_input_order='reversed',rejected_rehashed_false_exports=rejected,source_hashes={str(p):sha(p) for p in [Path(__file__),Path(producer.__file__),Path(reader.__file__)]},scientific_eligibility=False,scope=plan['scope'])
        a.output.parent.mkdir(exist_ok=True,parents=True)
        with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))


if __name__=='__main__':main()
