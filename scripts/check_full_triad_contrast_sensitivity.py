#!/usr/bin/env python3
"""Full-grid software checks, near-zero contracts and committed restart recovery."""
import argparse
import copy
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import tempfile
import summarize_full_triad_contrast_sensitivity as producer
import readback_full_triad_contrast_sensitivity as reader
from summarize_full_triad_order_robustness import summarize
from full_triad_robustness_sources import NUMERIC_FIELDS
from full_triad_contrast_sensitivity_sources import MASKS, CORES
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    with tempfile.TemporaryDirectory(prefix='fungal-full-contrast-contracts-') as temp:
        root=Path(temp);raw=root/'raw.tsv.gz';groups=root/'groups.jsonl.gz';triads=[];rows=[];group_records=[]
        screens=[dict(id='n50_c70'),dict(id='n30_c50')]
        for i in range(12):
            models=[[f'a{i}',1],[f'b{i}',2],[f'r{i}',3]]
            if i==1:models[0],models[1]=models[1],models[0]
            tid=hashlib.sha256(json.dumps(models,separators=(',',':')).encode()).hexdigest();triad=dict(triad_id=tid,models=models);triads.append(triad)
            for m,c in itertools.product(MASKS[:2],CORES[:2]):
                states=[]
                for order in itertools.product([0,1],repeat=3):
                    contrast=1.;status='computed_unique_at_numeric_tolerance'
                    if i==1:contrast=-1.
                    elif i==2:contrast=1. if m=='full' else -1.
                    elif i==3:contrast=1. if c=='reference_common' else -1.
                    elif i==4:contrast=-1. if order==(1,1,1) else 1.
                    elif i==5 and m=='plddt70' and c=='cycle_consistent' and order==(1,1,1):contrast=None;status='source_excluded'
                    elif i==6:contrast=.1;status='computed_nonunique_at_numeric_tolerance'
                    elif i==7:contrast=5e-10
                    elif i==8:contrast=1e-9
                    elif i==10:contrast=0.
                    elif i==11:contrast=1.001e-9
                    r=dict(triad_id=tid,mask=m,mapping_definition=c,fit_status=status,triples_sha256='synthetic-core',source_exclusions='')
                    r.update({field:1. for field in NUMERIC_FIELDS});r.update(rmsd_ar=2.+contrast if contrast is not None else '',rmsd_br=2. if contrast is not None else '',rmsd_ar_minus_br=contrast if contrast is not None else '')
                    for a,model in zip(['a','b','reference'],models):r['model_'+a],r['version_'+a]=model
                    for edge,o in zip(['ab','ar','br'],order):r['order_'+edge]=o
                    for s in screens:
                        sid=s['id'];core=int(status=='computed_unique_at_numeric_tolerance');inherited=int(not(i==9 and sid=='n50_c70'))
                        r.update({sid+'_pass':core*inherited,sid+'_core_pass':core,sid+'_three_pair_pass':inherited,sid+'_exclusions':'' if core*inherited else 'fixture_excluded'})
                    rows.append(r);states.append({k:str(v) for k,v in r.items()})
                group_records.append(dict(triad_id=tid,models=models,role_order=['a','b','reference'],mask=m,mapping_definition=c,
                    order_bits_layout=[list(o) for o in itertools.product([0,1],repeat=3)],**summarize(states,screens)))
        with gzip.open(raw,'wt') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(reversed(rows))
        with gzip.open(groups,'wt') as f:
            for g in group_records:f.write(json.dumps(g)+'\n')
        plan=dict(output=str(root/'production'),screens=screens,checkpoint_triads=2,contrast_numerical_tolerance_angstrom=1e-9,
            resources=dict(minimum_free_disk_gib=0,maximum_output_gib=1),
            expected=dict(measured_triads=12,source_fit_rows=384,source_order_groups=48,sensitivity_groups=108,screen_decisions=216,summary_rows=126),
            scope='Synthetic full-grid software contract; production source closure I/O stubbed, not biological evidence.')
        pp=root/'plan.json';pp.write_text(json.dumps(plan))
        def sources(config,path):return raw,groups,triads,{str(path):sha(path),str(raw):sha(raw),str(groups):sha(groups)}
        producer.load_sources=reader.load_sources=sources
        original=producer.scenario;calls=0
        def interrupted(*a,**kw):
            nonlocal calls
            calls+=1
            if calls==20:raise RuntimeError('owned_full_contrast_interruption')
            return original(*a,**kw)
        producer.scenario=interrupted
        try:producer.run(pp)
        except RuntimeError as e:assert str(e)=='owned_full_contrast_interruption'
        else:raise AssertionError('Missing owned interruption')
        finally:producer.scenario=original
        receipt=producer.run(pp);checked=reader.run(pp,root/'readback.json');assert receipt['checked_reused_chunks']==1
        assert checked['screen_decisions']==216
        with gzip.open(Path(plan['output'])/'contrast_sensitivity.jsonl.gz','rt') as f:exported=[json.loads(line) for line in f]
        joint=lambda i:exported[i*9+8]
        assert [joint(i)['numerical_tolerance_direction'] for i in range(12)]==['positive','negative','sign_uncertain','sign_uncertain','sign_uncertain','unavailable','nonunique_fit','within_numerical_tolerance','within_numerical_tolerance','positive','within_numerical_tolerance','positive']
        assert joint(8)['strict_zero_direction']=='positive' and joint(10)['strict_zero_direction']=='within_numerical_tolerance'
        assert joint(9)['screens']['n50_c70']['qualified_direction']=='excluded_by_quality' and joint(9)['screens']['n30_c50']['qualified_direction']=='positive'
        assert joint(5)['expected_orders']==32 and joint(5)['available_orders']==31
        mutations={
            'favorable_mask_reselection':lambda g,r:g[26].update(numerical_tolerance_direction='positive'),
            'favorable_core_reselection':lambda g,r:g[35].update(numerical_tolerance_direction='positive'),
            'dropped_unfavorable_order':lambda g,r:g[44].update(minimum_contrast_angstrom=1.),
            'promoted_incomplete_fit':lambda g,r:g[53].update(available_orders=32),
            'promoted_nonunique_fit':lambda g,r:g[62].update(all_selected_orders_unique=True),
            'nearzero_as_positive':lambda g,r:g[71].update(numerical_tolerance_direction='positive'),
            'changed_role_identity':lambda g,r:g[8]['models'].reverse(),
            'integer_unique_flag':lambda g,r:g[8].update(all_selected_orders_unique=1),
            'boolean_order_count':lambda g,r:g[8].update(available_orders=True),
            'invented_quality_pass':lambda g,r:g[89]['screens']['n50_c70'].update(qualified_direction='positive',all_selected_orders_pass=True),
            'changed_source_key':lambda g,r:g[8]['source_order_group_keys'][0].__setitem__(0,'wrong'),
            'missing_group':lambda g,r:g.pop(),
            'duplicate_group':lambda g,r:g.insert(0,copy.deepcopy(g[0])),
            'invented_scientific_eligibility':lambda g,r:r.update(scientific_eligibility=True),
            'changed_numerical_tolerance':lambda g,r:r.update(contrast_numerical_tolerance_angstrom=0.),
            'false_aggregate_count':lambda g,r:r['direction_counts'].update({'both_masks|both_cores|positive':999}),
        }
        rejected=[]
        for label,mutate in mutations.items():
            case=root/label;case.mkdir();config={**plan,'output':str(case)};cp=root/(label+'.json');cp.write_text(json.dumps(config))
            candidate=copy.deepcopy(exported);r=copy.deepcopy(receipt);mutate(candidate,r)
            with gzip.open(case/'contrast_sensitivity.jsonl.gz','wt') as f:
                for g in candidate:f.write(json.dumps(g)+'\n')
            chunks=case/'checkpoints';chunks.mkdir()
            for i,start in enumerate(range(0,len(candidate),18)):
                (chunks/f'{i:06d}.jsonl.gz').write_bytes(gzip.compress(('\n'.join(json.dumps(v) for v in candidate[start:start+18])+'\n').encode(),mtime=0))
            (case/'contrast_sensitivity_counts.tsv').write_bytes((Path(plan['output'])/'contrast_sensitivity_counts.tsv').read_bytes())
            bindings=sources(config,cp)[-1];configuration=case/'configuration.json';configuration.write_text(json.dumps(dict(plan_sha256=sha(cp),source_hashes=bindings)))
            r.update(plan_sha256=sha(cp),source_hashes={**bindings,str(configuration):sha(configuration)},artifacts={str(p.relative_to(case)):sha(p) for p in [case/'contrast_sensitivity.jsonl.gz',case/'contrast_sensitivity_counts.tsv',*chunks.glob('*.jsonl.gz')]})
            (case/'receipt.json').write_text(json.dumps(r))
            try:reader.run(cp,root/(label+'-readback.json'))
            except AssertionError:rejected.append(label)
            else:raise AssertionError('Accepted false full contrast export '+label)
        result=dict(status='passed_full_triad_contrast_sensitivity_software_contracts',**plan['expected'],checked_reused_chunks=1,
            raw_fit_input_order='reversed',rejected_rehashed_false_exports=rejected,source_hashes={str(p):sha(p) for p in [Path(__file__),Path(producer.__file__),Path(reader.__file__)]},
            scope=plan['scope'],scientific_eligibility=False)
        args.output.parent.mkdir(exist_ok=True,parents=True)
        with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))


if __name__=='__main__':main()
