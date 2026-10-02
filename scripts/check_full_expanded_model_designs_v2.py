#!/usr/bin/env python3
"""Software contracts for full recipe scope, exact deduplication and rank audits."""
import argparse
from collections import Counter
import copy
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from full_expanded_model_input_sources import schema, FIELDS, NUMERIC, INTEGER, ORDERS, MASKS, STRATUM, COUNT_FIELDS, DISTANCE, IDENTITY, NUISANCE
from prepare_full_expanded_model_designs import rank_audit, run
from readback_full_expanded_model_designs_v2 import numeric_audit
from run_ortholog_pair_guide_comparison import sha


def write(path,value):
    Path(path).parent.mkdir(parents=True,exist_ok=True);Path(path).write_text(json.dumps(value,indent=2)+'\n')


def table(path,rows,fields):
    with (gzip.open(path,'wt') if str(path).endswith('.gz') else Path(path).open('w')) as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def h(value):return hashlib.sha256(value.encode()).hexdigest()


def setup(root):
    for name in ['cases','covariance','selection','matching','inputs']: (root/name).mkdir()
    guides=['profile','mafft'];policies=['p1','p2'];scenarios=['S01','S54'];screens=[dict(id='screen1'),dict(id='screen2')]
    cases=[];rng=np.random.default_rng(82)
    for i in range(48):
        cases.append(dict(case_id=h('case'+str(47-i)),physical_case_id=h('physical'+str(i//2)),guide=guides[int(i>=24)],
            target_id=h('target'+str(i)),background_id=h('background'+str(i%3)),background_pair_key=h('background-pair'+str(i%2)),
            target_gene_a='distinct-gene-'+str(i),target_gene_b='gene-'+str(i+100),
            joint_full_pass_bits=3 if i%3 else 1,joint_plddt70_pass_bits=3 if i%4 else 1,
            joint_both_pass_bits=3 if i%3 and i%4 else 1))
    links=[]
    for policy,scenario in itertools.product(policies,scenarios):
        for i,c in enumerate(cases):
            if scenario=='S54' and i>=6:continue
            links.append(dict(source_row_ordinal=len(links)+1,case_id=c['case_id'],guide=c['guide'],physical_case_id=c['physical_case_id'],target_id=c['target_id'],background_id=c['background_id'],policy=policy,scenario_id=scenario))
    reuse=Counter(r['case_id'] for r in links)
    for c in cases:c['selection_records']=reuse[c['case_id']]
    table(root/'cases/case_index.tsv.gz',cases,list(cases[0]));table(root/'cases/selection_case_links.tsv.gz',links,list(links[0]))
    cov=[dict(case_id=c['case_id'],family_component=h('component'+str(i%2)),species_pattern_id=h('pattern'+str(i%5)),species_pattern_row=i%5) for i,c in enumerate(cases)]
    table(root/'covariance/case_covariance_index.tsv.gz',cov,list(cov[0]))
    write(root/'covariance-completion.json',dict(status='complete_verified_full_expanded_covariance',scientific_eligibility=False,trees=['tree1','tree2','tree3','tree4','tree5']))
    write(root/'covariance-plan.json',dict(output=str(root/'covariance'),trees=['tree1','tree2','tree3','tree4','tree5']))
    write(root/'selection/scenarios.json',[dict(scenario_id=s) for s in scenarios])
    write(root/'matching-plan.json',dict(selection=str(root/'selection')))
    write(root/'cases-plan.json',dict(output=str(root/'cases'),matching_plan=str(root/'matching-plan.json')))
    parts=[]
    for mask,order in itertools.product(MASKS,ORDERS):
        rows=[]
        for i,c in enumerate(cases):
            row={f:0. if f in NUMERIC else 0 if f in INTEGER else '' for f in FIELDS}
            row.update({k:c[k] for k in row if k in c});row.update(case_row=i,input_id=h('input'+c['case_id']+mask+order),mask=mask,order_contrast=order,numerical_usable=1,
                joint_mask_pass_bits=c['joint_'+mask+'_pass_bits'],joint_both_pass_bits=c['joint_both_pass_bits'])
            a,b=rng.uniform(.05,.95,2)
            for k,f in enumerate(DISTANCE,1):row[f]=a**k-b**k if not(mask=='plddt70' and c['guide']=='mafft') else 0.
            a,b=rng.uniform(.05,.95,2)
            for k,f in enumerate(IDENTITY,1):row[f]=a**k-b**k
            for f in NUISANCE:row[f]=float(rng.normal())
            if mask=='plddt70':row['aligned_plddt70_fraction_delta']=0.
            if c['guide']=='mafft':row['original_coverage_delta']=row['aligned_identity_delta']
            row['rmsd_delta']=float(rng.normal()) if i else 0.;row['native_tm_dissimilarity_delta']=.125 if c['guide']=='profile' else float(rng.normal())
            rows.append(row)
        folder=root/'inputs/inputs';folder.mkdir(exist_ok=True);fp=folder/(mask+'-'+order+'.parquet')
        pq.write_table(pa.Table.from_pylist(rows,schema=schema()),fp,compression='zstd')
        parts.append(dict(mask=mask,order_contrast=order,path=str(fp.relative_to(root/'inputs')),rows=len(cases),sha256=sha(fp)))
    write(root/'inputs/partition_manifest.json',parts);counts=[]
    for g,p,s,m,gate,screen,o in itertools.product(guides,policies,scenarios,MASKS,['mask','both_masks'],screens,ORDERS):
        selected=[r for r in links if (r['guide'],r['policy'],r['scenario_id'])==(g,p,s)];lookup={c['case_id']:c for c in cases};flag=1<<screens.index(screen)
        count=sum(bool(lookup[r['case_id']]['joint_'+(m if gate=='mask' else 'both')+'_pass_bits'] & flag) for r in selected)
        row={f:0 for f in COUNT_FIELDS};row.update(dict(zip(STRATUM,[g,p,s,m,gate,screen['id'],o])),quality_retained_numeric_records=count);counts.append(row)
    table(root/'inputs/setting_counts.tsv',counts,COUNT_FIELDS)
    inp=dict(output=str(root/'inputs'),cases_plan=str(root/'cases-plan.json'),covariance_plan=str(root/'covariance-plan.json'),covariance_completion=str(root/'covariance-completion.json'),
        expected=dict(scenarios=2,selected_records=len(links)),guides=guides,policies=policies,screens=screens)
    write(root/'inputs-plan.json',inp)
    artifacts={str(p.relative_to(root/'inputs')):sha(p) for p in (root/'inputs').rglob('*') if p.is_file()}
    summary=dict(logical_cases=48,selected_records=len(links),settings=len(counts),future_model_setting_rows=len(counts)*12)
    write(root/'inputs/receipt.json',dict(status='complete_full_expanded_model_inputs_pending_independent_readback',plan_sha256=sha(root/'inputs-plan.json'),scientific_eligibility=False,artifacts=artifacts,**summary))
    bindings={str(p):sha(p) for p in root.rglob('*') if p.is_file()}
    write(root/'inputs/completion_archive.json',dict(status='complete_verified_full_expanded_model_inputs_archive',services=[{},{}],summary=summary,source_hashes=bindings))
    write(root/'inputs-completion.json',dict(status='complete_verified_full_expanded_model_inputs',exact_process_journals_checked=2,bound_source_hashes=len(bindings),
        full_hash_archive=str(root/'inputs/completion_archive.json'),full_hash_archive_sha256=sha(root/'inputs/completion_archive.json'),producer_receipt=str(root/'inputs/receipt.json'),
        producer_receipt_sha256=sha(root/'inputs/receipt.json'),scientific_eligibility=False,**summary))
    plan=dict(inputs_completion=str(root/'inputs-completion.json'),inputs_plan=str(root/'inputs-plan.json'),expected=dict(logical_cases=48,selected_records=len(links),input_setting_rows=len(counts),model_setting_rows=len(counts)*12),
        trees=['tree1','tree2','tree3','tree4','tree5'],output=str(root/'output'),resources=dict(minimum_free_disk_gib=0),pins={},scope='Synthetic closed source and journals; software contracts only, not biological production acceptance.')
    write(root/'plan.json',plan);return root/'plan.json'


def invoke(plan, output, good=True):
    r=subprocess.run([sys.executable,'scripts/readback_full_expanded_model_designs_v2.py','--plan',str(plan),'--output',str(output)],capture_output=True,text=True)
    if good and r.returncode:raise RuntimeError(r.stdout+r.stderr)
    if not good:assert r.returncode!=0,'Altered design export accepted'


def mutate(root,name):
    out=root/'output';rp=out/'receipt.json';receipt=json.loads(rp.read_text())
    if name in ['cohort_order','cohort_alias_collapse']:
        manifest=json.loads((out/'cohort_manifest.json').read_text());entry=next(r for r in manifest if r['records']>2);fp=out/entry['path']
        with np.load(fp) as a:rows=a['case_rows'].copy()
        if name=='cohort_order':rows=rows[::-1]
        else:rows[1]=rows[0]
        np.savez_compressed(fp,case_rows=rows);entry['sha256']=sha(fp);write(out/'cohort_manifest.json',manifest);receipt['artifacts'][entry['path']]=sha(fp);receipt['artifacts']['cohort_manifest.json']=sha(out/'cohort_manifest.json')
    elif name in ['omit_s54','omit_empty_setting','setting_fit_changed','setting_degree_changed']:
        fp=out/'model_settings.tsv.gz'
        with gzip.open(fp,'rt') as f:r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames;rows=list(r)
        if name=='omit_s54':rows=[r for r in rows if r['scenario_id']!='S54']
        elif name=='omit_empty_setting':rows=[r for r in rows if r['disposition']!='empty_setting']
        elif name=='setting_fit_changed':rows[0]['fit_input_id']='wrong'
        else:rows[0]['degree']='4'
        table(fp,rows,fields);receipt['artifacts'][fp.name]=sha(fp)
    else:
        isfit=name in ['response_hash','constant_promoted','missing_fit','extra_fit','tree_changed','source_contract']
        fp=out/('unique_fit_inputs.jsonl' if isfit else 'unique_designs.jsonl');rows=[json.loads(l) for l in fp.open()]
        row=next((r for r in rows if r['disposition']=='constant_response_requires_review'),rows[0]) if name=='constant_promoted' else next(r for r in rows if r['records']>15)
        if name=='wrong_condition':
            row=next(r for r in rows if r['disposition']=='full_rank_design');row['normalized_condition_number']*=2
        elif name=='wrong_rank':row['rank']+=1
        elif name=='drop_nonzero_column':row['active_column_indices']=row['active_column_indices'][:-1]
        elif name=='wrong_scale':row['column_l2_after_maxabs'][0]*=2
        elif name=='loose_rank_tolerance':row['rank_tolerance']*=1e9
        elif name=='scaled_hash':row['scaled_design_sha256']='wrong'
        elif name=='raw_hash':row['raw_design_sha256']='wrong'
        elif name=='component_hash':row['family_components_sha256']='wrong'
        elif name=='physical_pair_hash':row['background_physical_pairs_sha256']='wrong'
        elif name=='ordered_input_hash':row['ordered_input_ids_sha256']='wrong'
        elif name=='wrong_predictor':row['predictor_columns'][1]='wrong'
        elif name=='response_hash':row['response_sha256']='wrong'
        elif name=='constant_promoted':row['disposition']='ready_for_working_covariance_fit'
        elif name=='missing_fit':rows.pop()
        elif name=='extra_fit':rows.append(copy.deepcopy(row))
        elif name=='tree_changed':row['trees'].pop()
        elif name=='source_contract':row['source_contract']='wrong'
        else:raise AssertionError(name)
        fp.write_text(''.join(json.dumps(r)+'\n' for r in rows));receipt['artifacts'][fp.name]=sha(fp)
    write(rp,receipt)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    changes=['cohort_order','cohort_alias_collapse','omit_s54','omit_empty_setting','setting_fit_changed','setting_degree_changed','wrong_rank','drop_nonzero_column','wrong_scale','loose_rank_tolerance','scaled_hash','raw_hash','component_hash','physical_pair_hash','ordered_input_hash','wrong_predictor','response_hash','constant_promoted','missing_fit','extra_fit','tree_changed','source_contract','wrong_condition']
    rng=np.random.default_rng(91);matrix=np.column_stack([np.ones(40),rng.normal(size=(40,7))]);matrix[:,2]*=1e-200;matrix[:,3]*=1e200
    audit=rank_audit(matrix,1);assert numeric_audit(matrix,1,audit)=='full_rank_design'
    rng_boundary=np.random.default_rng(17);x=rng_boundary.normal(size=40);z=rng_boundary.normal(size=40)
    boundary=np.column_stack([np.ones(40),x,x+1e-13*z]);boundary_audit=rank_audit(boundary,1)
    assert boundary_audit['disposition']=='rank_boundary_requires_review'
    assert numeric_audit(boundary,1,boundary_audit)=='rank_boundary_requires_review'
    with tempfile.TemporaryDirectory(prefix='full-design-contract-') as directory:
        root=Path(directory);plan=setup(root);result=run(plan);invoke(plan,root/'readback.json')
        assert result['model_setting_rows']==3840 and result['unique_fit_inputs']<result['model_setting_rows']
        assert result['fit_input_status_counts']['ready_for_working_covariance_fit']>0
        for status in ['empty_setting','rank_deficient_requires_review','sequence_axis_uninformative_requires_review','insufficient_residual_dimension','constant_response_requires_review']:assert result['fit_input_status_counts'][status]>0,status
        try:run(plan)
        except AssertionError:pass
        else:raise AssertionError('Completed restart accepted')
        snapshot={str(f.relative_to(root/'output')):f.read_bytes() for f in (root/'output').rglob('*') if f.is_file()}
        for change in changes:
            mutate(root,change);invoke(plan,root/'bad-readback.json',False)
            for f in (root/'output').rglob('*'):
                if f.is_file() and str(f.relative_to(root/'output')) not in snapshot:f.unlink()
            for name,raw in snapshot.items():(root/'output'/name).write_bytes(raw)
        root2=root/'interrupt';root2.mkdir();pp=setup(root2)
        try:run(pp,stop_after_cohorts=1)
        except InterruptedError:pass
        else:raise AssertionError('Interruption not exercised')
        run(pp);invoke(pp,root2/'readback.json')
    sources=[Path(__file__),*[Path('scripts')/n for n in ['full_expanded_model_design_sources.py','prepare_full_expanded_model_designs.py','readback_full_expanded_model_designs_v2.py']]]
    result=dict(status='passed_full_expanded_model_design_v2_software_contracts',rejected_rehashed_exports=changes,
        full_interrupt_replay_passed=True,completed_restart_refused=True,extreme_scale_qr_svd_check_passed=True,
        full_grid_retained_and_exact_cohorts_shared=True,source_and_journal_fixtures_synthetic=True,rank_boundary_condition_readback_passed=True,
        source_hashes={str(f):sha(f) for f in sources},scope='48 distinct gene-context cases; shared physical/model contexts, reused controls, 3,840 fixed model settings with exact reuse, empty/underspecified/rank-deficient/uninformative/constant-response dispositions. Full synthetic closed sources and journals check software only; not a pilot or production acceptance. Exact-zero nuisance columns remain explicit; no approximate term deletion or calibrated inference.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
