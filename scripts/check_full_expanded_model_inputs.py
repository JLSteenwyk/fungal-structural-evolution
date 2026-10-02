#!/usr/bin/env python3
"""Exercise full input contracts with synthetic prior handoffs, not a pilot."""
import argparse
from collections import Counter
import copy
import csv
import gzip
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import pyarrow as pa
import pyarrow.parquet as pq
from check_full_matched_coverage_cases import node
from index_full_matching_cases import metadata
from full_matching_case_sources import CASE_FIELDS
from full_expanded_model_input_sources import FIELDS,COUNT_FIELDS
from full_expanded_case_measurement_sources import FIELDS as CONTRAST_FIELDS
from join_full_expanded_case_measurements import joined
from run_ortholog_pair_guide_comparison import sha


def write(path,value):
    Path(path).parent.mkdir(exist_ok=True,parents=True);Path(path).write_text(json.dumps(value,indent=2)+'\n')


def table(path,rows,fields):
    Path(path).parent.mkdir(exist_ok=True,parents=True)
    with (gzip.open(path,'wt') if str(path).endswith('.gz') else Path(path).open('w')) as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def invoke(script,plan,output=None,success=True):
    command=[sys.executable,'scripts/'+script,'--plan',str(plan)]
    if output is not None:command+=['--output',str(output)]
    r=subprocess.run(command,capture_output=True,text=True)
    if success and r.returncode:raise RuntimeError(r.stdout+r.stderr)
    if not success:assert r.returncode!=0,'False input export accepted'


def setup(root):
    roots={x:root/x for x in ['cases','measurements','catalog','covariance','graph','matching','selection']}
    for p in roots.values():p.mkdir()
    targets={};backs={}
    for tid,g,a,b in [('T1','profile','A','B'),('T2','profile','C','D'),('T3','mafft','E','E'),('T4','mafft','F','G'),('T5','mafft','H','I'),('T6','profile','A','B'),('T7','mafft','J','K')]:targets[tid]=node(tid,g,a,b,True)
    for bid,g,a,b in [('B1','profile','L','M'),('B2','profile','N','O'),('B3','mafft','P','Q'),('B4','mafft','R','R'),('B5','mafft','S','U'),('B6','profile','L','M'),('B7','mafft','V','W')]:backs[bid]=node(bid,g,a,b)
    for records in [targets,backs]:
        for n in records.values():
            for e in ['a','b']:
                n['length_'+e]=100;n['mean_ca_plddt_'+e]=76.0 if e=='a' else 84.0;n['fraction_ca_plddt_below50_'+e]=0.1
            n['sequence_distance']=0 if n['same_model'] else 1.75 if n['node_id'].startswith('T') else 0.25
    catalog={}
    for role,nodes in [('target',targets),('background',backs)]:
        for n in nodes.values():
            if n['same_model']:continue
            for mask,order in itertools.product(['full','plddt70'],[0,1]):
                key=role,n['pair_key'],mask,order
                if key in catalog:continue
                good=not(n['node_id']=='T5' and mask=='plddt70' and order==1)
                size=20 if n['node_id']=='T7' else 60+order*20
                value=dict(role=role,pair_key=n['pair_key'],mask=mask,order=order,numerical_usable=int(good),native_status='aligned',numerical_exclusion_reasons='' if good else 'quarantined_fixture',
                    model_a=n['model_id_a'],version_a=n['version_a'],model_b=n['model_id_b'],version_b=n['version_b'],original_length_a=100,original_length_b=100,
                    aligned_length=size,original_coverage_a=size/100,original_coverage_b=size/100,rmsd_recomputed=0 if role=='target' and n['node_id']=='T2' else 2+order if role=='target' else 3+order,
                    tm_left_native=.7-order*.1,tm_right_native=.6-order*.1,sequence_identity_exact=(.2+.6*order) if role=='target' else (.3+.2*order),joint_plddt70_fraction=1 if mask=='plddt70' else .5+order*.2,
                    pair_mask_pass_bits=0,both_masks_pass_bits=0)
                catalog[key]=value
    for role,nodes in [('target',targets),('background',backs)]:
        for n in nodes.values():
            for mask in ['full','plddt70']:
                good=not n['same_model'] and all(catalog[role,n['pair_key'],mask,o]['numerical_usable'] and catalog[role,n['pair_key'],mask,o]['aligned_length']>=30 for o in [0,1])
                n[mask+'_bits']=int(good)
            n['both_bits']=n['full_bits'] & n['plddt70_bits']
            if not n['same_model']:
                for mask,o in itertools.product(['full','plddt70'],[0,1]):
                    catalog[role,n['pair_key'],mask,o]['pair_mask_pass_bits']=n[mask+'_bits'];catalog[role,n['pair_key'],mask,o]['both_masks_pass_bits']=n['both_bits']
    for role,nodes in [('target',targets),('background',backs)]:
        (roots['graph']/(role+'_nodes.jsonl')).write_text(''.join(json.dumps(n)+'\n' for n in nodes.values()))
    pairs=[('T1','B1'),('T2','B1'),('T1','B2'),('T3','B3'),('T4','B4'),('T5','B5'),('T6','B6'),('T7','B7')]
    cases=[];selections=[];screens=[dict(id='n30_c50',minimum_aligned_residues=30,minimum_original_coverage=.5),dict(id='n50_c70',minimum_aligned_residues=50,minimum_original_coverage=.7)]
    for tid,bid in pairs:
        t,b=targets[tid],backs[bid];c=metadata(t,b)
        for mask in ['full','plddt70','both']:
            c['target_'+mask+'_pass_bits']=t[mask+'_bits'];c['control_'+mask+'_pass_bits']=b[mask+'_bits'];c['joint_'+mask+'_pass_bits']=t[mask+'_bits'] & b[mask+'_bits']
        scenarios=['S01'] if tid=='T1' and bid=='B1' else ['S54'] if tid=='T1' else ['S01','S54']
        for policy,scenario in itertools.product(['p1','p2'],scenarios):
            selections.append(dict(source_row_ordinal=len(selections)+1,case_id=c['case_id'],physical_case_id=c['physical_case_id'],target_id=tid,background_id=bid,guide=t['guide'],policy=policy,scenario_id=scenario))
        c.update(selection_records=2*len(scenarios),endpoint_order_bits=3,policy_bits=3,scenario_bits=3);cases.append(c)
    table(roots['cases']/'case_index.tsv.gz',cases,CASE_FIELDS)
    table(roots['cases']/'selection_case_links.tsv.gz',selections,list(selections[0]))
    catalog_strings={key:{k:str(v) for k,v in value.items()} for key,value in catalog.items()}
    table(roots['catalog']/'directed_measurements.tsv.gz',list(catalog.values()),list(next(iter(catalog.values()))))
    cs=[{k:str(v) for k,v in c.items()} for c in cases]
    contrast=[joined(c,m,catalog_strings) for c in cs for m in ['full','plddt70']]
    table(roots['measurements']/'case_mask_contrasts.tsv.gz',contrast,CONTRAST_FIELDS)
    covariance=[]
    for i,c in enumerate(cases):covariance.append({**{k:c[k] for k in ['case_id','physical_case_id','target_id','background_id','guide','target_family','background_family','selection_records']},'family_component':'shared-component','species_pattern_id':'pattern-'+c['guide'],'species_pattern_row':int(c['guide']=='mafft')})
    table(roots['covariance']/'case_covariance_index.tsv.gz',covariance,list(covariance[0]))
    scenarios=[dict(scenario_id='S01'),dict(scenario_id='S54')];write(roots['selection']/'scenarios.json',scenarios)
    attrs=[]
    for g,p,s,m in itertools.product(['profile','mafft'],['p1','p2'],['S01','S54'],['full','plddt70','both']):
        members=[r for r in selections if (r['guide'],r['policy'],r['scenario_id'])==(g,p,s)];ci={c['case_id']:c for c in cases}
        for bit,screen in enumerate(screens):attrs.append(dict(guide=g,policy=p,scenario_id=s,mask=m,screen=screen['id'],all_target_records=20,matched_records=len(members),unmatched_records=20-len(members),joint_pass_matched_records=sum(bool(ci[r['case_id']]['joint_'+m+'_pass_bits'] & (1<<bit)) for r in members)))
    table(roots['matching']/'matched_attrition_counts.tsv',attrs,list(attrs[0]))
    mp=root/'matching-plan.json';write(mp,dict(output=str(roots['matching']),graph=str(roots['graph']),selection=str(roots['selection'])))
    definitions={'cases':'complete_verified_full_matching_logical_case_index','measurements':'complete_verified_full_expanded_matched_case_measurements','catalog':'complete_verified_full_expanded_directed_measurement_catalog','covariance':'complete_verified_full_expanded_covariance'}
    paths={k:root/(k+'-plan.json') for k in definitions};completions={k:root/(k+'-completion.json') for k in definitions}
    for k in definitions:
        config=dict(output=str(roots[k]))
        if k in ['cases','catalog']:config['matching_plan']=str(mp)
        if k=='cases':config.update(guides=['profile','mafft'],policies=['p1','p2'],screens=screens)
        if k=='measurements':config.update(case_index_completion=str(completions['cases']),catalog_completion=str(completions['catalog']))
        if k=='covariance':config['case_completion']=str(completions['cases'])
        write(paths[k],config)
    count=dict(logical_cases=8,physical_cases=len({c['physical_case_id'] for c in cases}),selected_records=len(selections),unmatched_decisions=100,directed_states=len(catalog),scenarios=2)
    raw=[p for p in root.rglob('*') if p.is_file()]
    bindings={str(p):sha(p) for p in raw}
    for k,status in definitions.items():
        outputs=[p for p in roots[k].iterdir() if p.is_file()];rp=roots[k]/'receipt.json'
        write(rp,dict(status='synthetic_prior_receipt',plan_sha256=sha(paths[k]),scientific_eligibility=False,artifacts={p.name:sha(p) for p in outputs}))
        stage_bind={**bindings,str(rp):sha(rp)};summary={x:count[x] for x in ['logical_cases','physical_cases','selected_records']} if k!='catalog' else dict(directed_states=count['directed_states'])
        if k=='measurements':summary['unmatched_decisions']=count['unmatched_decisions']
        ap=root/(k+'-archive.json');write(ap,dict(status=status+'_archive',services=[{},{}],summary=summary,source_hashes=stage_bind))
        write(completions[k],dict(status=status,exact_process_journals_checked=2,bound_source_hashes=len(stage_bind),full_hash_archive=str(ap),full_hash_archive_sha256=sha(ap),producer_receipt=str(rp),producer_receipt_sha256=sha(rp),scientific_eligibility=False,**summary))
    plan=dict(expected=count,output=str(root/'output'),guides=['profile','mafft'],policies=['p1','p2'],screens=screens,pins={},resources=dict(minimum_free_disk_gib=0),scope='Synthetic prior source/journal handoffs; exercises full software contracts only, not scientific production acceptance.')
    for k in definitions:plan[k+'_plan']=str(paths[k]);plan[k+'_completion']=str(completions[k])
    pp=root/'plan.json';write(pp,plan);return pp


def mutate(root,name):
    pp=root/'output'/'partition_manifest.json';parts=json.loads(pp.read_text());part=next(p for p in parts if p['mask']=='full' and p['order_contrast']=='mean');path=root/'output'/part['path']
    rows=pq.read_table(path).to_pylist();row=rows[0]
    if name=='nonlinear_mean_after_transform':row['aligned_identity_square_delta']=row['aligned_identity_delta']**2
    elif name=='gene_distance_sign':row['gene_distance_delta']=-row['gene_distance_delta']
    elif name=='response_shift':row['rmsd_delta']+=1
    elif name=='confidence_shift':row['mean_ca_plddt_delta']+=1
    elif name=='coverage_shift':row['original_coverage_delta']+=.1
    elif name=='gene_alias_collapsed':row['target_gene_a']=rows[1]['target_gene_a']
    elif name=='species_pattern_changed':row['species_pattern_row']+=1
    elif name=='component_changed':row['family_component']='wrong'
    elif name=='state_key_swapped':row['target_state_keys']=row['background_state_keys']
    elif name=='input_identity_changed':row['input_id']='wrong'
    elif name=='source_case_row_changed':row['case_row']+=1
    elif name=='null_as_zero':rows[3]['rmsd_delta']=0.
    elif name=='invalid_promoted':rows[3]['numerical_usable']=1
    elif name=='drop_case':rows.pop()
    elif name=='extra_case':rows.append(copy.deepcopy(row))
    elif name in ['wrong_setting_count','omit_final_scenario','reuse_weight_changed']:
        cp=root/'output'/'setting_counts.tsv'
        with cp.open() as f:counts=list(csv.DictReader(f,delimiter='\t'))
        if name=='wrong_setting_count':counts[0]['quality_retained_records']=str(int(counts[0]['quality_retained_records'])+1)
        elif name=='reuse_weight_changed':counts[0]['sum_inverse_background_node_reuse_weights']='999'
        else:counts=[c for c in counts if c['scenario_id']!='S54']
        table(cp,counts,COUNT_FIELDS);receipt=json.loads((root/'output'/'receipt.json').read_text());receipt['artifacts'][cp.name]=sha(cp);write(root/'output'/'receipt.json',receipt);return
    else:raise AssertionError(name)
    pq.write_table(pa.Table.from_pylist(rows,schema=pq.read_schema(path)),path,compression='zstd');part['sha256']=sha(path);write(pp,parts)
    rp=root/'output'/'receipt.json';receipt=json.loads(rp.read_text());receipt['artifacts'][part['path']]=sha(path);receipt['artifacts'][pp.name]=sha(pp);write(rp,receipt)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    changes=['nonlinear_mean_after_transform','gene_distance_sign','response_shift','confidence_shift','coverage_shift','gene_alias_collapsed','species_pattern_changed','component_changed','state_key_swapped','input_identity_changed','source_case_row_changed','null_as_zero','invalid_promoted','drop_case','extra_case','wrong_setting_count','omit_final_scenario','reuse_weight_changed']
    with tempfile.TemporaryDirectory(prefix='expanded-model-input-contract-') as directory:
        base=Path(directory)/'base';base.mkdir();plan=setup(base);invoke('prepare_full_expanded_model_inputs.py',plan);invoke('readback_full_expanded_model_inputs.py',plan,base/'readback.json')
        invoke('prepare_full_expanded_model_inputs.py',plan,success=False)
        for change in changes:
            snapshot={str(p.relative_to(base)):p.read_bytes() for p in (base/'output').rglob('*') if p.is_file()}
            mutate(base,change);invoke('readback_full_expanded_model_inputs.py',plan,base/'bad-readback.json',success=False)
            for p in (base/'output').rglob('*'):
                if p.is_file() and str(p.relative_to(base)) not in snapshot:p.unlink()
            for p,raw in snapshot.items():(base/p).write_bytes(raw)
        replay=Path(directory)/'replay';replay.mkdir();rp=setup(replay)
        command=[sys.executable,'-c',"import sys;sys.path.insert(0,'scripts');from prepare_full_expanded_model_inputs import run;run(sys.argv[1],stop_after_cases=1)",str(rp)]
        attempt=subprocess.run(command,capture_output=True,text=True);assert attempt.returncode and 'Software interruption contract' in attempt.stderr
        invoke('prepare_full_expanded_model_inputs.py',rp);invoke('readback_full_expanded_model_inputs.py',rp,replay/'readback.json')
    result=dict(status='passed_full_expanded_model_input_software_contracts',rejected_rehashed_exports=changes,
        full_interrupt_replay_passed=True,completed_restart_refused=True,source_and_journal_fixtures_synthetic=True,
        source_hashes={str(p):sha(p) for p in [Path(__file__),'scripts/full_expanded_model_input_sources.py','scripts/prepare_full_expanded_model_inputs.py','scripts/readback_full_expanded_model_inputs.py']},
        scope='Eight gene-context cases including shared model/physical aliases, retained true zero, target/control same-model exclusions, a one-order quarantine, coverage exclusions, original S54, background reuse and both missing/zero distinctions. Complete fixture setting grid and independent Decimal/SQLite checks. Native case/catalog/covariance/measurement closure proofs and prior journals are synthetic software contracts, not a biological pilot or production acceptance.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
