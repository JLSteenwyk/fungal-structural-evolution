#!/usr/bin/env python3
"""Complete synthetic-grid contracts for uniform covariance qualification."""
import argparse
from collections import Counter
import csv
import gzip
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from check_full_entity_operators import setup as bank_setup
from prepare_full_entity_operators import run as bank_run
from readback_full_entity_operators import run as bank_readback
from full_entity_operator_sources import SUMMARY_FIELDS as BANK_SUMMARY
from full_expanded_model_input_sources import FIELDS,NUMERIC,INTEGER,schema,ORDERS,MASKS,NUISANCE
from full_expanded_model_design_sources import AXES,OUTCOMES,DEGREES,SETTING_FIELDS,array_digest,cohort_id,design_id,fit_id,digest
from prepare_full_expanded_model_designs import rank_audit
from prepare_full_covariance_qualification import run
from run_ortholog_pair_guide_comparison import sha
from scipy import sparse
from covariance_basis_context import ComponentKernelProducts
from covariance_basis_independent import IndependentKernelProducts
from full_covariance_qualification_sources import covariance_status
from readback_full_covariance_qualification import numeric


def write(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')


def close_fixture(root,status,producer,summary,additional):
    p=json.loads(Path(producer).read_text());bindings=dict(p.get('source_hashes',{}));bindings.update(additional)
    bindings[str(producer)]=sha(producer)
    for name,d in p['artifacts'].items():bindings[str(Path(producer).parent/name)]=d
    archive=root/'completion_archive.json'
    write(archive,dict(status=status+'_archive',source_hashes=bindings,services=[{},{}],summary=summary))
    completion=root/'completion.json'
    write(completion,dict(status=status,scientific_eligibility=False,full_hash_archive=str(archive),
        full_hash_archive_sha256=sha(archive),bound_source_hashes=len(bindings),exact_process_journals_checked=2,
        producer_receipt=str(producer),producer_receipt_sha256=sha(producer),**summary))
    return completion


def setup(root):
    bankdir=root/'bank';bankdir.mkdir();bankplan=bank_setup(bankdir);b=bank_run(bankplan)
    bank_readback(bankplan,bankdir/'output/readback.json')
    bank_summary={k:b[k] for k in BANK_SUMMARY}
    bankcompletion=close_fixture(bankdir,'complete_verified_full_entity_operator_bank',bankdir/'output/receipt.json',bank_summary,
        {str(bankdir/'output/readback.json'):sha(bankdir/'output/readback.json'),str(bankplan):sha(bankplan)})
    ids=json.loads((bankdir/'output/case_ids.json').read_text());n=len(ids);trees=b['trees']
    inputs=root/'inputs';inputs.mkdir();rng=np.random.default_rng(582);parts=[];arrays={}
    for mask,order in itertools.product(MASKS,ORDERS):
        records=[]
        for i,cid in enumerate(ids):
            row={k:0. if k in NUMERIC else 0 if k in INTEGER else '' for k in FIELDS}
            row.update(case_id=cid,input_id=digest(['fixture-input',cid,mask,order]),case_row=i,numerical_usable=1,
                mask=mask,order_contrast=order)
            for axis,columns in AXES.items():
                a,b=rng.uniform(.05,.95,2)
                for degree,column in enumerate(columns,1):row[column]=float(a**degree-b**degree)
            for column in NUISANCE:row[column]=float(rng.normal())
            if mask=='plddt70':row['aligned_plddt70_fraction_delta']=0.
            row['rmsd_delta']=float(rng.normal());row['native_tm_dissimilarity_delta']=.125 if mask=='plddt70' else float(rng.normal())
            records.append(row)
        fp=inputs/(mask+'-'+order+'.parquet');pq.write_table(pa.Table.from_pylist(records,schema=schema()),fp)
        parts.append(dict(mask=mask,order_contrast=order,path=fp.name,rows=n,sha256=sha(fp)))
        arrays[mask,order]=records
    write(inputs/'partition_manifest.json',parts)
    inputplan=root/'inputs-plan.json';write(inputplan,dict(output=str(inputs)))
    input_summary=dict(logical_cases=n,model_input_rows=n*10)
    write(inputs/'receipt.json',dict(status='complete_full_expanded_model_inputs_pending_independent_readback',
        source_hashes={str(inputplan):sha(inputplan)},artifacts={p.name:sha(p) for p in inputs.iterdir() if p.is_file()},**input_summary))
    inputcompletion=close_fixture(inputs,'complete_verified_full_expanded_model_inputs',inputs/'receipt.json',input_summary,{})
    designs=root/'designs';designs.mkdir();folder=designs/'cohorts';folder.mkdir();cohorts=[]
    for number,mask,indices in [(1,'full',range(n)),(2,'full',range(12)),(3,'plddt70',range(12,n)),
                                (4,'full',range(6)),(5,'full',[]),(6,'plddt70',range(2))]:
        rows=np.asarray(sorted(indices,key=lambda i:ids[i]),dtype=np.int64)
        ordered=array_digest([ids[i] for i in rows],'S64');cid=cohort_id('profile',mask,ordered)
        fp=folder/(cid+'.npz');np.savez_compressed(fp,case_rows=rows)
        cohorts.append(dict(cohort_id=cid,guide='profile',mask=mask,ordered_case_ids_sha256=ordered,
            case_rows_sha256=array_digest(rows,'<i8'),path=str(fp.relative_to(designs)),sha256=sha(fp),
            records=len(rows),membership_occurrences=2,fixture_scenario='S'+str(number),members=rows))
    cohorts.sort(key=lambda c:c['cohort_id']);settings=[];design_records=[];fit_records=[];ds=Counter();fs=Counter()
    contract=digest(['synthetic-input-design-contract'])
    for cohort in cohorts:
        rows=cohort['members'];cid=cohort['cohort_id']
        for order,axis,degree in itertools.product(ORDERS,AXES,DEGREES):
            values=arrays[cohort['mask'],order];columns=['intercept',*AXES[axis][:degree],*NUISANCE]
            matrix=np.column_stack([np.ones(len(rows)),*[np.asarray([values[i][k] for i in rows]) for k in columns[1:]]])
            audit=rank_audit(matrix,degree);did=design_id(cid,order,axis,degree,contract)
            design=dict(design_id=did,cohort_id=cid,mask=cohort['mask'],order_contrast=order,sequence_axis=axis,degree=degree,
                predictor_columns=columns,raw_design_sha256=array_digest(matrix,'<f8'),
                ordered_input_ids_sha256=array_digest([values[i]['input_id'] for i in rows],'S64'),**audit)
            design_records.append(design);ds[design['disposition']]+=1
            for outcome in OUTCOMES:
                y=np.asarray([values[i][outcome] for i in rows]);fid=fit_id(did,outcome,array_digest(y,'<f8'))
                status=design['disposition'] if design['disposition']!='full_rank_design' else 'constant_response_requires_review' if np.ptp(y)==0 else 'ready_for_working_covariance_fit'
                fit_records.append(dict(fit_input_id=fid,design_id=did,outcome=outcome,disposition=status));fs[status]+=1
                for policy in ['p1','p2']:
                    settings.append(dict(guide='profile',policy=policy,scenario_id=cohort['fixture_scenario'],mask=cohort['mask'],
                        eligibility_gate='mask',screen='screen1',order_contrast=order,outcome=outcome,sequence_axis=axis,degree=degree,
                        cohort_id=cid,design_id=did,fit_input_id=fid,records=len(rows),disposition=status,nominal_tree_fits=5))
    write(designs/'cohort_manifest.json',[{k:v for k,v in c.items() if k not in ['members','fixture_scenario']} for c in cohorts])
    for name,records in [('unique_designs.jsonl',design_records),('unique_fit_inputs.jsonl',fit_records)]:
        (designs/name).write_text(''.join(json.dumps(r)+'\n' for r in records))
    with gzip.open(designs/'model_settings.tsv.gz','wt') as f:
        w=csv.DictWriter(f,SETTING_FIELDS,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(settings)
    designplan=root/'design-plan.json';write(designplan,dict(output=str(designs),trees=trees))
    design_summary=dict(logical_cases=n,model_setting_rows=len(settings),unique_cohorts=len(cohorts),unique_designs=len(design_records),trees=trees)
    write(designs/'receipt.json',dict(status='complete_full_expanded_model_designs_pending_independent_readback',
        source_hashes={str(designplan):sha(designplan)},artifacts={str(p.relative_to(designs)):sha(p) for p in designs.rglob('*') if p.is_file()},**design_summary))
    designcompletion=close_fixture(designs,'complete_verified_full_expanded_model_designs',designs/'receipt.json',design_summary,{})
    plan=root/'plan.json';write(plan,dict(design_completion=str(designcompletion),operator_completion=str(bankcompletion),inputs_completion=str(inputcompletion),
        design_plan=str(designplan),operator_plan=str(bankplan),inputs_plan=str(inputplan),trees=trees,expected=dict(logical_cases=n,model_setting_rows=len(settings)),
        output=str(root/'output'),pins={},resources=dict(minimum_free_disk_gib=0),scope='Synthetic source closures and journals; complete software contracts only, no biological pilot or production acceptance.'))
    return plan


def invoke(plan,output,good=True):
    result=subprocess.run([sys.executable,'scripts/readback_full_covariance_qualification.py','--plan',str(plan),'--output',str(output)],capture_output=True,text=True)
    if good and result.returncode:raise RuntimeError(result.stdout+result.stderr)
    if not good:assert result.returncode!=0,'Altered covariance qualification accepted'


def mutate(root,name):
    output=root/'output';receipt=json.loads((output/'receipt.json').read_text())
    if name.startswith('link_'):
        fp=output/'setting_audit_links.tsv.gz'
        with gzip.open(fp,'rt') as f:r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames;records=list(r)
        if name=='link_omit_final_tree':records=[r for r in records if r['tree']!='t5']
        elif name=='link_omit_empty':records=[r for r in records if r['disposition']!='empty_setting']
        elif name=='link_wrong_fit':records[0]['fit_input_id']='wrong'
        elif name=='link_wrong_audit':records[0]['audit_id']='wrong'
        elif name=='link_wrong_status':records[0]['combined_disposition']='qualified_uniform_working_covariance_basis'
        else:raise AssertionError(name)
        with gzip.open(fp,'wt') as f:w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(records)
    else:
        fp=output/'design_covariance_audits.jsonl.gz'
        with gzip.open(fp,'rt') as f:records=[json.loads(line) for line in f]
        row=next(r for r in records if r['numerical_audit'] is not None);value=row['numerical_audit']
        if name=='omit_final_tree':records=[r for r in records if r['tree']!='t5']
        elif name=='omit_unsigned':records=[r for r in records if r['loading_mode']!='unsigned']
        elif name=='wrong_raw_gram':value['raw_gram'][0][0]*=2
        elif name=='wrong_projected_gram':value['projected_gram'][0][0]*=2
        elif name=='wrong_envelope':value['projected_roundoff_envelope'][0][0]*=2
        elif name=='wrong_rank':value['reml_diagnostics']['rank']+=1
        elif name=='promote_dependent_basis':row['disposition']='qualified_uniform_working_covariance_basis'
        elif name=='undo_target_fold':row['folded_terms']['target_node']='separately_estimated'
        elif name=='wrong_family_fold':row['folded_terms']['family']='independent_endpoint_and_intercept_variances'
        elif name=='wrong_design_hash':row['raw_design_sha256']='wrong'
        elif name=='wrong_source_contract':row['source_contract']='wrong'
        elif name=='false_failure':row.update(disposition='numerical_covariance_qualification_requires_review',numerical_audit=None,error_type='ArithmeticError',error_message='fake')
        else:raise AssertionError(name)
        with gzip.open(fp,'wt') as f:f.write(''.join(json.dumps(r)+'\n' for r in records))
    receipt['artifacts'][fp.name]=sha(fp);write(output/'receipt.json',receipt)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    # Exercise the qualified path as well as the dependent production fixture.
    # Random block-local bases are a numerical contract, not biological data.
    rng=np.random.default_rng(194);labels=np.repeat(np.arange(4),6);operators={}
    for name in ['background_node','model_pair','gene','model']:
        z=np.zeros((24,12))
        for block in range(4):z[block*6:(block+1)*6,block*3:(block+1)*3]=rng.normal(size=(6,3))
        operators[name]=sparse.csr_matrix(z)
    operators['family_intercept']=sparse.csr_matrix((np.ones(24),(np.arange(24),labels)),shape=(24,4))
    factor=rng.normal(size=(24,4));x=np.column_stack([np.ones(24),rng.normal(size=(24,2))])
    value=ComponentKernelProducts(labels,operators,np.ones(24)).tree(factor).audit(x)
    assert covariance_status(value)=='qualified_uniform_working_covariance_basis'
    numeric(value,IndependentKernelProducts(labels,operators).tree(factor).project(x),value['kernel_names'],24)
    changes=['omit_final_tree','omit_unsigned','wrong_raw_gram','wrong_projected_gram','wrong_envelope','wrong_rank','promote_dependent_basis',
        'undo_target_fold','wrong_family_fold','wrong_design_hash','wrong_source_contract','false_failure','link_omit_final_tree',
        'link_omit_empty','link_wrong_fit','link_wrong_audit','link_wrong_status']
    with tempfile.TemporaryDirectory(prefix='full-uniform-covariance-contract-') as directory:
        root=Path(directory);plan=setup(root);result=run(plan);invoke(plan,root/'readback.json')
        assert result['model_setting_rows']==720 and result['audit_rows']==1800 and result['setting_audit_links']==7200
        for status in ['empty_setting','insufficient_residual_dimension','covariance_basis_requires_review']:assert result['audit_status_counts'][status]>0
        assert result['link_status_counts']['constant_response_requires_review']>0
        try:run(plan)
        except AssertionError:pass
        else:raise AssertionError('Completed qualification restart accepted')
        snapshot={str(p.relative_to(root/'output')):p.read_bytes() for p in (root/'output').rglob('*') if p.is_file()}
        for change in changes:
            mutate(root,change);invoke(plan,root/'bad-readback.json',False)
            for name,raw in snapshot.items():(root/'output'/name).write_bytes(raw)
        root2=root/'interrupted';root2.mkdir();pp=setup(root2)
        try:run(pp,stop_after_cohorts=1)
        except InterruptedError:pass
        else:raise AssertionError('Interruption not exercised')
        run(pp);invoke(pp,root2/'readback.json')
    sources=[Path(__file__),*[Path('scripts')/v for v in ['full_covariance_qualification_sources.py','prepare_full_covariance_qualification.py',
        'readback_full_covariance_qualification.py','covariance_basis_context.py','covariance_basis_independent.py','covariance_basis_audit.py']]]
    result=dict(status='passed_full_uniform_covariance_qualification_contracts',synthetic_model_setting_rows=720,
        synthetic_design_audits=1800,synthetic_setting_links=7200,rejected_rehashed_exports=changes,
        full_interrupt_replay_passed=True,completed_restart_refused=True,source_and_journal_fixtures_synthetic=True,
        independent_seven_basis_qualified_path_passed=True,
        all_loading_modes_trees_empty_and_review_cases_retained=True,source_hashes={str(p):sha(p) for p in sources},
        scope='Complete 24-case/six-cohort synthetic closed-grid software contracts, with signed/unsigned five-tree '
              'products, exact target/family folds, independent latent-space/SQL replay and all original settings. '
              'No production qualification, biological pilot, variance fit or calibrated inference.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
