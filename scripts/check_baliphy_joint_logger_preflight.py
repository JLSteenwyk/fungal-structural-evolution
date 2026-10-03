#!/usr/bin/env python3
"""Full future startup grid, three capped native fixtures and serialized contracts."""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from unittest.mock import patch

from ancestral_chain_attempt import run_attempt, sha, write_json
from prepare_baliphy_joint_logger_preflight import build
import run_baliphy_joint_logger_preflight as producer
import readback_baliphy_joint_logger_preflight as reader


def rejected(action):
    try:action()
    except (AssertionError,KeyError,ValueError,FileNotFoundError):return
    raise AssertionError('Invalid startup accepted')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);a=p.parse_args();root=a.output.resolve();root.mkdir(exist_ok=False)
    future_path=Path('metadata/baliphy_joint_node_logger_future_models_20261003_v4.json')
    future=json.loads(future_path.read_text());producer.verify(future)
    native=json.loads(Path(future['source_plan']).read_text());producer.verify(native)
    roles=json.loads(Path(future['future_roles']).read_text());originals=json.loads(Path(native['jobs']).read_text())
    binary=Path(originals[0]['config']['command'][5]);prlimit=Path('/usr/bin/prlimit')
    jobs=build(roles,originals,binary,prlimit,[]);assert len(jobs)==1620
    alterations=[]
    for case in ['missing_role','duplicate_id','duplicate_source_id','duplicate_seed','reused_seed',
                 'invalid_seed','changed_role','changed_alias','changed_prior','changed_program_hash','claimed_execution']:
        changed=copy.deepcopy(roles)
        if case=='missing_role':changed.pop()
        elif case=='duplicate_id':changed[1]['chain']['chain_id']=changed[0]['chain']['chain_id']
        elif case=='duplicate_source_id':changed[1]['source_chain_id']=changed[0]['source_chain_id']
        elif case=='duplicate_seed':changed[1]['chain']['seed']=changed[0]['chain']['seed']
        elif case=='reused_seed':changed[0]['chain']['seed']=originals[0]['chain']['seed']
        elif case=='invalid_seed':changed[0]['chain']['seed']=0
        elif case=='changed_role':changed[0]['chain']['chain']=4
        elif case=='changed_alias':changed[0]['chain']['original_configuration_ids']=['invented']
        elif case=='changed_prior':changed[0]['chain']['prior_label']='invented'
        elif case=='changed_program_hash':changed[0]['chain']['program_sha256']='0'*64
        else:changed[0]['native_execution_launched']=True
        rejected(lambda:build(changed,originals,binary,prlimit,[]));alterations.append(case)
    fixture=Path('data/software_audits/baliphy-joint-node-logger-20261003-v4').resolve()
    bindings={str(future_path):sha(future_path),str(Path(future['future_roles'])):sha(future['future_roles']),
              str(Path(native['jobs'])):sha(native['jobs'])}
    actual={};configs={};receipts={}
    for index,prior in enumerate(['broad','centered','package']):
        old=fixture/(prior+'-original_tip_logger.hs');new=fixture/(prior+'-joint_node_logger.hs')
        alignment=fixture/'alignment.faa';tree=fixture/'tree.nwk'
        config=dict(command=[str(prlimit),'--as='+str(8*2**30),'--cpu=120','--fsize='+str(64*2**20),
            '--',str(binary),'--seed',str(20267001+index),'run',str(new),'--test','--log-format','json'],
            timeout_seconds=180,pins={str(path):sha(path) for path in [prlimit,binary,old,new,alignment,tree]})
        receipt=run_attempt(root/'native'/prior,config)
        chain=dict(chain_id='software-'+prior,effective_input_group='software',chain=1,prior_label=prior,
            family='synthetic',proteins=5,original_configuration_ids=['software-'+prior],seed=20267001+index,
            program=str(old),program_sha256=sha(old),alignment=str(alignment),alignment_sha256=sha(alignment),
            tree=str(tree),tree_sha256=sha(tree))
        job=dict(chain=chain,fresh_seed=20267001+index,program=str(new),config=config)
        row=producer.inspect(job,receipt,'software')
        assert row['status']=='reference_startup_homology_density_and_representation_checked',row
        actual[prior]=row;configs[prior]=config;receipts[prior]=receipt
        for path in [old,new,alignment,tree,receipt,receipt.parent.parent/'configuration.json']:bindings[str(path)]=sha(path)
        for name,h in json.loads(receipt.read_text())['artifacts'].items():bindings[str(receipt.parent/name)]=h
        print('actual_joint_logger_native_startup',prior,row['status'],flush=True)
    # Full metadata and actual source/config generation above; only the following
    # full producer/reader run mocks native execution and per-role source admission.
    mocked_jobs=copy.deepcopy(jobs);fake={}
    for index,job in enumerate(mocked_jobs):
        c=job['chain'];prior=c['prior_label'];job['config']=configs[prior]
        row=copy.deepcopy(actual[prior]);row.update(chain_id=c['chain_id'],effective_input_group=c['effective_input_group'],
            model_input_identity=c['effective_input_group']+'-'+prior,chain_role=c['chain'],family=c['family'],
            original_configuration_ids=c['original_configuration_ids'],source_seed=c['seed'],fresh_seed=job['fresh_seed'])
        if index in [0,4]:row.update(status='unsuccessful_native_startup_retained',exit_code=1)
        fake[c['chain_id']]=row
    stage=root/'artificial-full-stage';job_path=root/'artificial-jobs.json';write_json(job_path,mocked_jobs)
    plan_path=root/'artificial-plan.json';plan=dict(jobs=str(job_path),output=str(stage),pins={},
        resources=dict(memory_gib=16,minimum_free_disk_gib=1),scope='Explicitly mocked native/source-admission software fixture only')
    write_json(plan_path,plan)
    def fake_inspect(job,receipt,digest):
        row=copy.deepcopy(fake[job['chain']['chain_id']]);row['plan_sha256']=digest;return row
    serial=[]
    with patch.object(producer,'run_attempt',side_effect=lambda root,config:receipts[next(k for k,v in configs.items() if v==config)]),\
         patch.object(producer,'inspect',side_effect=fake_inspect),patch.object(reader,'inspect',side_effect=fake_inspect),\
         patch.object(sys,'argv',['software','--plan',str(plan_path)]):
        producer.main();result=json.loads((stage/'receipt.json').read_text())
        assert result['validated_startups']==1618 and result['unresolved_startup_quartets']==2
        rejected(producer.main);serial.append('producer_restart')
        for case in ['missing_checkpoint','changed_seed','scientific_acceptance','false_posterior']:
            cp=next((stage/'chains').glob('*.json'));original=cp.read_bytes()
            if case=='missing_checkpoint':cp.unlink()
            else:
                row=json.loads(original)
                if case=='changed_seed':row['fresh_seed']+=1
                elif case=='scientific_acceptance':row['scientific_eligibility']=True
                else:row['posterior_sampling_launched']=True
                write_json(cp,row)
            rejected(lambda:reader.readback(plan_path));cp.write_bytes(original);serial.append(case)
        reread=reader.readback(plan_path);assert reread['validated_startups']==1618
        rejected(lambda:reader.readback(plan_path));serial.append('reader_restart')
    for name in ['prepare_baliphy_joint_logger_preflight','run_baliphy_joint_logger_preflight',
        'readback_baliphy_joint_logger_preflight','check_baliphy_joint_logger_preflight',
        'baliphy_joint_node_logger_v3','baliphy_full_node_logger','baliphy_reference_initialization',
        'baliphy_reference_startup_readback_v2','run_baliphy_reference_preflight']:
        path=Path('scripts',name+'.py');bindings[str(path)]=sha(path)
    receipt=dict(status='passed_full_joint_logger_preflight_software_contracts',
        checked_utc=datetime.now(timezone.utc).isoformat(),full_roles=1620,full_quartets=405,
        full_source_grid_checked=True,native_startups=3,native_fixture_priors=list(actual),
        actual_native_startup_checks=list(actual.values()),full_producer_reader_serialization_checked=True,
        artificial_failures_retained=2,artificial_unresolved_quartets_retained=2,
        altered_designs_rejected=alterations,serialization_rejections=serial,source_hashes=bindings,
        scientific_eligibility=False,posterior_sampling_launched=False,
        scope='All1620future roles and405programs/source inverse/configuration generation checked. Three capped native synthetic--test initializations validate reference homology, fixed-tip/free-ancestor representation, degree-aware densities and runtime clades across allpriors. Full1620producer/readerserialization uses explicitly mocked native execution/source admission andretains2failures/2unresolvedquartets. No fungal pilot,MCMC,actualfullgridstartup,jointframes,memoryrepair,posteriorqualification or scientificacceptance.')
    with a.receipt.open('x') as h:h.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_hashes','actual_native_startup_checks']}),flush=True)


if __name__=='__main__':main()
