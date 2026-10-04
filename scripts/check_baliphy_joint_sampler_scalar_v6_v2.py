#!/usr/bin/env python3
"""Qualify full V6 role construction and real retained output admission.

Reuses three original V6 native fixtures, never restarts them. Altered records
and a synthetic native failure exist only in fresh private fixture copies.
Full-grid summary accounting is explicitly synthetic, not 1620 native runs.
"""
import argparse
import ast
import copy
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

from ancestral_chain_attempt import sha
from baliphy_joint_sampler_scalar_v6 import build_jobs,inspect,summarize,SUCCESS,INVALID,REVIEW
from check_baliphy_scalar_json_logger_v6_v9 import seed_census
from reference_measurement_union_sources import bind,verify


def rejected(action):
    try:action()
    except (AssertionError,ValueError,KeyError,TypeError):return
    raise AssertionError('Altered role construction accepted')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists();root=a.output.resolve();root.mkdir(exist_ok=False)
    future_path=Path('metadata/baliphy_scalar_json_v6_future_models_20261004_v1.json')
    future=json.loads(future_path.read_text());verify(future['pins'])
    assert future['status']=='prepared_full_future_scalar_json_v6_models_not_launched'
    old_path=Path('metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json')
    old=json.loads(old_path.read_text());verify(old['pins'])
    roles=json.loads(Path(future['future_roles']).read_text());previous=json.loads(Path(old['jobs']).read_text())
    pins=dict(future['pins']);pins.update(old['pins'])
    # Bind the complete transitive project source imports before construction
    # and rehash them after every actual/altered-output check. V1 verified the
    # old plan up front but did not flatten these inherited provider bindings.
    pending=[Path(__file__),Path('scripts/baliphy_joint_sampler_scalar_v6.py')]
    dependencies=set()
    while pending:
        module=pending.pop().resolve()
        if module in dependencies:continue
        dependencies.add(module);bind(pins,module)
        for node in ast.walk(ast.parse(module.read_text())):
            names=[node.module] if isinstance(node,ast.ImportFrom) and node.module else [a.name for a in node.names] if isinstance(node,ast.Import) else []
            for name in names:
                local=Path('scripts')/(name.split('.')[0]+'.py')
                if local.exists():pending.append(local)
    forbidden=seed_census(pins)
    assert len(forbidden)==future['forbidden_seed_count']
    jobs=build_jobs(roles,previous,forbidden)
    job_path=root/'software_grid_jobs.json';job_path.write_text(json.dumps(jobs,indent=2)+'\n')
    design_negatives=[]
    for case in ['missing_role','duplicate_id','duplicate_source','duplicate_seed','old_seed','bool_seed',
                 'changed_prior','changed_alias','changed_role','wrong_schema_seed_namespace',
                 'changed_program_digest','invented_native_execution']:
        changed=copy.deepcopy(roles)
        if case=='missing_role':changed.pop()
        elif case=='duplicate_id':changed[1]['chain']['chain_id']=changed[0]['chain']['chain_id']
        elif case=='duplicate_source':changed[1]['source_v5_chain_id']=changed[0]['source_v5_chain_id']
        elif case=='duplicate_seed':changed[1]['chain']['seed']=changed[0]['chain']['seed']
        elif case=='old_seed':changed[0]['chain']['seed']=previous[0]['chain']['seed']
        elif case=='bool_seed':changed[0]['chain']['seed']=True
        elif case=='changed_prior':changed[0]['chain']['prior_label']='invented'
        elif case=='changed_alias':changed[0]['chain']['original_configuration_ids']=[]
        elif case=='changed_role':changed[0]['chain']['chain']=5
        elif case=='wrong_schema_seed_namespace':changed[0]['seed_namespace']='old'
        elif case=='changed_program_digest':changed[0]['chain']['program_sha256']='0'*64
        else:changed[0]['native_execution_launched']=True
        rejected(lambda:build_jobs(changed,previous,forbidden));design_negatives.append(case)
    gate_path=Path('metadata/baliphy_scalar_json_v6_software_validation_20261004_v9.json')
    gate=json.loads(gate_path.read_text());verify(gate['source_hashes'])
    mapping=root/'fixture_mapping.tsv'
    with mapping.open('x') as f:
        writer=csv.DictWriter(f,delimiter='\t',fieldnames=['guide','family','dataset','level','source_node','retained_set_json'])
        writer.writeheader()
        for i,n in enumerate([2,3,4,5]):writer.writerow(dict(guide='profile',family='software',dataset='whole',level=i,source_node='n'+str(i),retained_set_json=json.dumps(list('abcde')[:n])))
    fixture=Path('data/software_audits/baliphy-joint-node-logger-20261003-v4').resolve()
    actual={};fixture_jobs={};fixture_receipts={}
    for item in gate['paired_prior_checks']:
        prior=item['prior'];directory=Path(item['new_directory']);receipt=directory.parent/'receipt.json'
        config_path=receipt.parent.parent/'configuration.json';config=json.loads(config_path.read_text())
        program=config['command'][config['command'].index('run')+1]
        base=next(j['chain'] for j in jobs if j['chain']['prior_label']==prior)
        chain=dict(base,chain_id='software-v6-'+prior,seed=item['seed'],family='software',proteins=5,
                   effective_input_group='software',original_configuration_ids=['software-whole-'+prior],
                   alignment=str(fixture/'alignment.faa'),alignment_sha256=sha(fixture/'alignment.faa'),
                   tree=str(fixture/'tree.nwk'),tree_sha256=sha(fixture/'tree.nwk'),program=program,program_sha256=sha(program))
        job=dict(chain=chain,config=config,source_seed=20267001,memory_reservation_bytes=12*2**30)
        export=root/'actual_fixture_arrays'/prior
        row=inspect(job,receipt,'software-v6',mapping,export)
        assert row['status']==SUCCESS and row['scalar_integrity_accepted']
        assert row['scalar_v6_audit']['rows']==21 and row['scalar_v6_audit']['mapped_values_compared']==903
        assert len(row['joint_frames'])==3
        assert inspect(job,receipt,'software-v6',mapping,export,False)==row
        # A readback with absent exports must fail, without recreating them.
        missing=root/'absent_export'/prior
        rejected(lambda:inspect(job,receipt,'software-v6',mapping,missing,False));assert not missing.exists()
        actual[prior]=row;fixture_jobs[prior]=job;fixture_receipts[prior]=receipt
        bind(pins,receipt);bind(pins,config_path)
        for path,h in config['pins'].items():bind(pins,path,h)
        for path,h in json.loads(receipt.read_text())['artifacts'].items():bind(pins,receipt.parent/path,h)
    # Every change is confined to a new copy of one actual native attempt.
    # Rehash output receipts so scalar rejection is exercised beyond custody.
    cases=[];private_rows={}
    for case in ['missing_json','old_header','duplicate_iteration','duplicate_key','changed_finite_value',
                 'wrong_quality_count','bare_infinity','wrong_placeholder',
                 'tagged_infinity_review','literal_null_review','synthetic_native_failure']:
        original=fixture_receipts['broad'];folder=root/'altered'/case
        shutil.copytree(original.parent.parent,folder)
        receipt=folder/'attempt-0001'/'receipt.json';directory=next(receipt.parent.glob('independent-chain-*'))
        path=directory/'C1.log.json';lines=path.read_text().splitlines()
        header=json.loads(lines[0]);row=json.loads(lines[1]);parameters=row['parameters//']['S1/'];field='ASRV.Gamma:alpha'
        assert field in parameters and type(parameters[field]) in (int,float)
        quality=row['numericParameterQuality//']
        if case=='missing_json':path.unlink()
        else:
            if case=='old_header':header.pop('projectScalarSchema');lines[0]=json.dumps(header)
            elif case=='duplicate_iteration':row['iter']=1
            elif case=='changed_finite_value':parameters[field]*=10
            elif case=='wrong_quality_count':quality['numericLeafCount']+=1
            elif case=='wrong_placeholder':
                parameters[field]='__project_scalar_v6__:negative_infinity'
                quality['nonfinite']=[dict(path=['S1/',field],kind='positive_infinity')]
            elif case=='tagged_infinity_review':
                parameters[field]='__project_scalar_v6__:positive_infinity'
                quality['nonfinite']=[dict(path=['S1/',field],kind='positive_infinity')]
            elif case=='literal_null_review':
                parameters[field]='__project_scalar_v6__:literal_null';quality['literalNullPaths']=[['S1/',field]]
                quality['numericLeafCount']-=1
            lines[1]=json.dumps(row)
            if case=='duplicate_key':lines[1]=lines[1].replace('"iter": 0','"iter": 0, "iter": 0',1)
            if case=='bare_infinity':lines[1]=lines[1].replace('"'+field+'": '+str(parameters[field]),'"'+field+'": Infinity',1)
            path.write_text('\n'.join(lines)+'\n')
        native=json.loads(receipt.read_text())
        native['artifacts']={str(p.relative_to(receipt.parent)):sha(p) for p in receipt.parent.rglob('*') if p.is_file() and p!=receipt}
        if case=='synthetic_native_failure':native.update(exit_code=1,status='failed')
        receipt.write_text(json.dumps(native,indent=2)+'\n')
        export=root/'altered_exports'/case
        row=inspect(fixture_jobs['broad'],receipt,'software-v6',mapping,export)
        expected=REVIEW if case in ['tagged_infinity_review','literal_null_review'] else 'unsuccessful_sampler_qualification_attempt_retained' if case=='synthetic_native_failure' else INVALID
        assert row['status']==expected and not row['scalar_integrity_accepted'] and not row['joint_frames']
        assert not export.exists()
        cases.append(dict(case=case,status=row['status'],native_failure_is_synthetic=case=='synthetic_native_failure'))
        private_rows[case]=row
    # Synthetic whole-grid accounting of retained actual fixture templates.
    # No role above is substituted for a real fungal execution or output.
    synthetic=[]
    for i,job in enumerate(jobs):
        chain=job['chain'];template=private_rows['changed_finite_value'] if i==0 else private_rows['tagged_infinity_review'] if i==4 else actual[chain['prior_label']]
        row=copy.deepcopy(template)
        row.update(chain_id=chain['chain_id'],effective_input_group=chain['effective_input_group'],
                   model_input_identity=chain['effective_input_group']+'-'+chain['prior_label'],
                   chain_role=chain['chain'],prior_label=chain['prior_label'],family=chain['family'],
                   original_configuration_ids=chain['original_configuration_ids'],seed=chain['seed'],source_seed=job['source_seed'])
        synthetic.append(row)
    summary=summarize(synthetic)
    assert summary['full_chains']==1620 and summary['complete_quartets']==403 and summary['unresolved_quartets']==2
    assert summary['scalar_review_roles']==1 and summary['scalar_v6_finite_checked_roles']==1618
    assert summary['joint_saved_frames']==4854 and summary['scalar_v6_mapped_values_checked']==1618*903
    for path in [Path(__file__),Path('scripts/baliphy_joint_sampler_scalar_v6.py'),future_path,old_path,
                 Path(old['jobs']),Path(future['future_roles']),gate_path,mapping,job_path]:bind(pins,path)
    for path in root.rglob('*'):
        if path.is_file():bind(pins,path)
    verify(pins)
    result=dict(status='passed_v6_scalar_joint_admission_adapter_software_v2',checked_utc=datetime.now(timezone.utc).isoformat(),
        full_source_roles=1620,full_source_models=405,effective_inputs=135,original_configuration_aliases=324,
        forbidden_seed_count=len(forbidden),transitive_project_source_modules=len(dependencies),
        inherited_v3_plan_provider_pins_merged=True,design_cases_rejected=design_negatives,
        retained_actual_native_roles=3,retained_actual_scalar_rows=63,retained_actual_mapped_values=2709,
        retained_actual_joint_frames=9,readback_missing_export_cases_rejected=3,
        rehashed_private_output_cases=cases,synthetic_full_grid_summary=summary,
        source_hashes=pins,new_native_runs=0,full_grid_native_execution_launched=False,
        full_workflow_serialization_qualified=False,full_grid_v6_startup_qualified=False,
        scientific_eligibility=False,posterior_qualified=False,gpu=False,new_cost_usd=0,
        scope='Complete original1620role construction plus three retained actual V6 native output admissions '
              '(63scalar rows,2709mapped values,9joint frames). Twelve altered role constructions rejected; '
              'eleven rehashed private output cases retain invalid/nonfinite/null/native-failure states '
              'without arrays. Full1620role summary accounting is synthetic fixture-template expansion, '
              'not actual executions or full controller/closure qualification. No original output, '
              'installed software, prior/model, production horizon or job changed.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','synthetic_full_grid_summary']},indent=2))


if __name__=='__main__':main()
