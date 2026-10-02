#!/usr/bin/env python3
"""Check full-grid selection, recovery, coefficient provenance and false-export rejection."""
import argparse
from collections import Counter
from contextlib import redirect_stdout
import copy
from datetime import datetime,timezone
import gzip
import io
from itertools import combinations
import json
from pathlib import Path
import tempfile
from unittest.mock import patch
import export_full_whole_protein_comparisons as producer
import readback_full_whole_protein_comparisons as reader
from run_whole_protein_ml import digest
from screen_duplication_alignment_reuse import sha

PASS='ml_candidate_passed_numerical_optimization_checks'
FLAG='ml_candidate_requires_optimization_review'
ERROR='fit_error_requires_review'


def write(path,value):
    Path(path).parent.mkdir(exist_ok=True,parents=True);Path(path).write_text(json.dumps(value,allow_nan=False)+'\n')


def fixture(root):
    fitplan=root/'fit_plan.json';write(fitplan,dict(scope='software fixture only'))
    trees=['tree_'+str(i) for i in range(5)];original={};follow={};ids=[];specs=[];support_rows={};bindings={str(fitplan):sha(fitplan)}
    follow_root=root/'followup';follow_root.mkdir()
    for i in range(10):
        cols=['outcome','intercept','sequence']
        if i in [1,2,6]:cols.append('length')
        if i==7:cols[-1]='coverage'
        spec=dict(columns=cols,records=90 if i==9 else 100,fixture_case=i)
        identifier=digest(spec);ids.append(identifier);specs.append(spec);width=len(cols)-1
        scales=[2.0]*(width-1);beta=[float(j+1) for j in range(width)];cov=[[float(j==k)*4 for k in range(width)] for j in range(width)]
        conversion=[1.]+[.5]*(width-1)
        params=dict(beta=beta,conditional_beta_covariance=cov,raw_unit_beta=[beta[j]*conversion[j] for j in range(width)],
            raw_unit_conditional_beta_covariance=[[cov[j][k]*conversion[j]*conversion[k] for k in range(width)] for j in range(width)],profiled_scale=2.0,variance_profile_denominator=spec['records'])
        state=FLAG if i in [2,3,4,5] else ERROR if i==6 else PASS
        ih=digest(dict(input=i));support_rows[identifier]=dict(input_sha256=ih,geometry_id='geometry_'+str(i),classification='numerical_support_review' if i==7 else 'zero_supported_to_numeric_tolerance',original_classification='original_review',source_kind='fixture_support')
        for t in trees:
            payload=dict(status=state,error='retained') if state==ERROR else dict(status=state,negative_profiled_ml=[100.,90.,91.,96.,105.,98.,0.,97.,100.,50.][i],covariate_scales=scales,log1p_ratios=[.1,.2,.3],**params)
            saved=dict(fit_input_id=identifier,tree=t,input_sha256=ih,specification=spec,payload=payload,payload_sha256=digest(payload),plan_sha256=sha(fitplan))
            op=root/'original'/identifier/(t+'.json');write(op,saved)
            entry=dict(fit_input_id=identifier,tree=t,path=str(op),sha256=sha(op),status=state);original[identifier,t]=entry;bindings[str(op)]=sha(op)
            if state==FLAG:
                fs='numerical_followup_requires_review' if i==4 else 'numerical_followup_error_requires_review' if i==5 else 'numerical_followup_passed_pending_independent_readback'
                result=dict(status=fs,error='retained numerical failure') if i==5 else dict(status=fs,fitted=dict(**params,negative_profiled_ml=89. if i==2 else 95. if i==3 else 104.),selected_theta=[.2,.3,.4],recovery=None if i==2 else dict(selection='fixture'),raw_unit_beta=params['raw_unit_beta'],raw_unit_conditional_beta_covariance=params['raw_unit_conditional_beta_covariance'])
                case=dict(identity=dict(fit_input_id=identifier,tree=t,original_fit=str(op),original_fit_sha256=sha(op),input_sha256=ih),result=result,result_sha256=digest(result),original_status=FLAG,scientific_eligibility=False)
                cp=follow_root/identifier/(t+'.json');write(cp,case)
                follow[identifier,t]=dict(fit_input_id=identifier,tree=t,path=str(cp.relative_to(follow_root)),sha256=sha(cp),status=fs);bindings[str(cp)]=sha(cp)
    class Support:
        def resolve(self,identifier,ih):
            r=support_rows[identifier];assert r['input_sha256']==ih;return r
    links=[]
    for a,b in combinations(range(10),2):
        left,right=set(specs[a]['columns'][1:]),set(specs[b]['columns'][1:])
        relation='different_observations_no_direct_comparison' if specs[a]['records']!=specs[b]['records'] else 'identical_named_design' if left==right else 'left_named_columns_nested_in_right' if left<right else 'right_named_columns_nested_in_left' if right<left else 'same_observations_no_named_column_nesting'
        links.append(dict(comparison_id=str(a)+'_'+str(b),left_input=ids[a],right_input=ids[b],relation=relation))
    link_root=root/'links';link_root.mkdir();lp=link_root/'comparison_input_map.jsonl';lp.write_text(''.join(json.dumps(r)+'\n' for r in links));bindings[str(lp)]=sha(lp)
    lr=dict(relation_counts=dict(Counter(r['relation'] for r in links)))
    plan=dict(pins={},fit_plan=str(fitplan),links=str(link_root),output=str(root/'export'),expected=dict(unique_inputs=10,fit_summaries=50,parameter_records=50,trees=5,comparisons_per_tree=45,comparisons=225),resources=dict(minimum_free_disk_gib=0,maximum_output_gib=1),scope='Software full five-tree/selection grid only; production full75070 inputs375350 fits4147200 comparisons and closure journals are not stubbed or claimed complete.')
    pp=root/'plan.json';write(pp,plan);bindings[str(pp)]=sha(pp)
    def sources(plan,path):return copy.deepcopy(original),copy.deepcopy(follow),follow_root,lr,Support(),dict(bindings)
    return pp,sources,ids


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();rejected=[]
    with tempfile.TemporaryDirectory(prefix='full-whole-protein-comparison-') as td,redirect_stdout(io.StringIO()):
        root=Path(td);plan,sources,ids=fixture(root);out=root/'export'
        class Interrupted(Exception):pass
        real_usage=producer.shutil.disk_usage;calls=0
        def stop_after_first_checkpoint(path):
            nonlocal calls
            calls+=1
            if calls==2:raise Interrupted()
            return real_usage(path)
        with patch.object(producer,'load_sources',side_effect=sources),patch.object(reader,'load_sources',side_effect=sources):
            with patch.object(producer.shutil,'disk_usage',side_effect=stop_after_first_checkpoint):
                try:producer.run(plan)
                except Interrupted:pass
                else:raise AssertionError('Checkpoint interruption was not exercised')
            assert len(list(out.glob('*.checkpoint.json')))==1 and not (out/'receipt.json').exists()
            parameter_hash=sha(out/'selected_fit_parameters.jsonl.gz');r=producer.run(plan)
            assert r['checked_reused_tree_checkpoints']==1 and sha(out/'selected_fit_parameters.jsonl.gz')==parameter_hash
            verified=reader.run(plan,out/'readback.json')
            try:producer.run(plan)
            except AssertionError:completed_restart_refused=True
            else:raise AssertionError('Completed export restarted')
            summaries=[json.loads(l) for l in (out/'fit_summaries.jsonl').read_text().splitlines()]
            by_id={r['fit_input_id']:r for r in summaries if r['tree']=='tree_0'}
            assert by_id[ids[2]]['negative_profiled_ml']==89. and by_id[ids[2]]['selected_source_kind']=='full_grid_refinement'
            assert by_id[ids[3]]['negative_profiled_ml']==95. and by_id[ids[3]]['selected_source_kind']=='full_grid_recovery'
            assert by_id[ids[4]]['effective_status']=='ml_full_grid_followup_requires_review'
            assert by_id[ids[5]]['effective_status']==FLAG and by_id[ids[5]]['selected_source_kind']=='original_retained_after_followup_error'
            assert by_id[ids[6]]['negative_profiled_ml'] is None and by_id[ids[6]]['numerical_fit_verified'] is False
            originals,follow,follow_root,lr,support,bindings=sources({},plan)
            entry=originals[ids[0],'tree_0'];saved=json.loads(Path(entry['path']).read_text())
            for label,change in [
                ('source_coefficient_unit_transform',lambda s:s['payload'].update(raw_unit_beta=[99.,99.])),
                ('source_covariance_unit_transform',lambda s:s['payload'].update(raw_unit_conditional_beta_covariance=[[99.,0.],[0.,99.]])),
                ('source_nonpositive_covariate_scale',lambda s:s['payload'].update(covariate_scales=[0.])),
                ('source_wrong_parameter_dimension',lambda s:s['payload'].update(beta=[1.]))]:
                changed=copy.deepcopy(saved);change(changed);changed['payload_sha256']=digest(changed['payload'])
                try:reader.reconstruct(entry,changed,follow,follow_root,support)
                except (AssertionError,ValueError):rejected.append(label)
                else:raise AssertionError('Invalid raw source parameter accepted: '+label)
            baseline={p.name:p.read_bytes() for p in out.iterdir() if p.is_file()}
            def lines(name):
                raw=gzip.decompress((out/name).read_bytes()).decode() if name.endswith('.gz') else (out/name).read_text()
                return [json.loads(l) for l in raw.splitlines()]
            def save_lines(name,rows):
                raw=''.join(json.dumps(v,allow_nan=False)+'\n' for v in rows).encode();(out/name).write_bytes(gzip.compress(raw,mtime=0) if name.endswith('.gz') else raw)
            def corrupt(name,mutate,label):
                for n,data in baseline.items():(out/n).write_bytes(data)
                rows=lines(name);mutate(rows);save_lines(name,rows)
                # Rehash all output evidence, so rejection must come from source/arithmetic contracts.
                config=dict(plan_sha256=sha(plan),summary_sha256=sha(out/'fit_summaries.jsonl'),parameters_sha256=sha(out/'selected_fit_parameters.jsonl.gz'))
                for checkpoint in out.glob('*.checkpoint.json'):
                    c=json.loads(checkpoint.read_text());c['configuration']=config;c['sha256']=sha(out/(checkpoint.stem.replace('.checkpoint','')+'.jsonl.gz'));write(checkpoint,c)
                rr=json.loads((out/'receipt.json').read_text());rr['artifacts']={n:sha(out/n) for n in rr['artifacts']};write(out/'receipt.json',rr)
                try:reader.run(plan,out/(label+'.json'))
                except (AssertionError,KeyError,ValueError,TypeError):rejected.append(label)
                else:raise AssertionError('False export accepted: '+label)
            sf='fit_summaries.jsonl';pf='selected_fit_parameters.jsonl.gz';tf='tree_0.jsonl.gz'
            def item(rows,i):return next(r for r in rows if r['fit_input_id']==ids[i] and r['tree']=='tree_0')
            corrupt(sf,lambda rows:item(rows,4).update(effective_status='ml_refinement_passed_numerical_checks'),'promoted_unresolved_followup')
            corrupt(sf,lambda rows:item(rows,5).update(effective_status='ml_refinement_passed_numerical_checks'),'promoted_followup_error')
            corrupt(sf,lambda rows:item(rows,6).update(negative_profiled_ml=1.,numerical_fit_verified=True),'promoted_original_error')
            corrupt(sf,lambda rows:item(rows,3).update(selected_source_path='false_source'),'wrong_selected_source')
            corrupt(sf,lambda rows:item(rows,2).update(followup_source_sha256='0'*64),'wrong_followup_hash')
            corrupt(sf,lambda rows:item(rows,2).update(negative_profiled_ml=91.),'retained_unselected_original_objective')
            corrupt(sf,lambda rows:item(rows,7).update(support_classification='zero_supported_to_numeric_tolerance'),'promoted_support_review')
            corrupt(sf,lambda rows:item(rows,6).update(numerical_fit_verified=0),'boolean_integer_equivalence')
            corrupt(sf,lambda rows:rows.pop(),'missing_original_summary')
            corrupt(pf,lambda rows:item(rows,3).update(raw_unit_beta=[9.,9.]),'changed_raw_coefficients')
            corrupt(pf,lambda rows:item(rows,2).update(scaled_conditional_beta_covariance=[[9.]*3]*3),'changed_conditional_covariance')
            corrupt(pf,lambda rows:item(rows,2).update(covariate_scales=[3.,3.]),'changed_covariate_scale')
            corrupt(pf,lambda rows:item(rows,3).update(selected_log1p_ratios=[0.,0.,0.]),'changed_selected_variance_ratios')
            corrupt(pf,lambda rows:item(rows,6).update(raw_unit_beta=[1.,1.,1.]),'manufactured_error_coefficients')
            corrupt(pf,lambda rows:item(rows,2).update(scientific_eligibility=True),'promoted_scientific_eligibility')
            corrupt(pf,lambda rows:rows.pop(),'missing_parameter_record')
            corrupt(tf,lambda rows:rows[0].update(log_likelihood_gain_right_vs_left=0.),'changed_likelihood_gain')
            corrupt(tf,lambda rows:rows[0].update(nested_added_fixed_coefficients=1.0),'integer_float_equivalence')
            corrupt(tf,lambda rows:rows[0].update(nesting_status='false_nesting'),'changed_nesting_status')
            corrupt(tf,lambda rows:rows[0].update(comparison_status='scientifically_accepted'),'promoted_comparison')
            corrupt(tf,lambda rows:rows.pop(),'missing_comparison_row')
            for n,data in baseline.items():(out/n).write_bytes(data)
            rr=json.loads((out/'receipt.json').read_text());rr['followup_cases_integrated']=5;write(out/'receipt.json',rr)
            try:reader.run(plan,out/'false_followup_census.json')
            except AssertionError:rejected.append('false_followup_census')
            else:raise AssertionError('False full follow-up census accepted')
        proof=dict(software_inputs=10,software_trees=5,software_fit_summaries=50,software_parameter_records=50,software_followups=20,software_comparisons=225,
            relation_counts=r['relation_counts'],effective_fit_status_counts=r['effective_fit_status_counts'],selected_source_kind_counts=r['selected_source_kind_counts'],
            checked_reused_tree_checkpoints=1,deterministic_parameter_restart_hash=True,completed_restart_refused=completed_restart_refused,rejected_false_exports=len(rejected),rejected_cases=rejected)
    sources_paths=[Path(__file__),Path(producer.__file__),Path(reader.__file__),Path('scripts/full_whole_protein_comparison_sources.py'),Path('scripts/export_whole_protein_ml_comparisons.py'),Path('scripts/readback_whole_protein_ml_comparisons.py')]
    result=dict(status='passed_full_whole_protein_comparison_software_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),**proof,
        source_hashes={str(p.relative_to(Path.cwd()) if p.is_absolute() else p):sha(p) for p in sources_paths},scientific_eligibility=False,
        scope='Software contracts only. Full five-tree fixture and every source status/refinement/recovery/error branch use real producer/independent selection and comparison reader. Closed source/journal I/O replaced only for fixture. Production full75070 inputs375350 fits4147200 comparisons never stubbed or claimed complete. Coefficient covariance conditional, no inferential acceptance.')
    a.output.parent.mkdir(exist_ok=True,parents=True)
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
