#!/usr/bin/env python3
"""Full software workflow with mock execution and retained native fixture replay."""
import argparse
import copy
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import resource
import time
from unittest.mock import patch

import numpy as np

from ancestral_chain_attempt import sha, write_json
from baliphy_joint_sampler_qualification import build_jobs, inspect, summarize, SUCCESS, INVALID
from baliphy_joint_sampler_gates import prerequisites, resource_boundaries
from independent_joint_ancestral_frames import verify_arrays, write_arrays
from independent_short_sampler_outputs_v2 import strict_json
from run_baliphy_reference_preflight import verify
import run_baliphy_joint_sampler_qualification as workflow


def rejected(action):
    try: action()
    except (AssertionError,ValueError,KeyError,RuntimeError,FileNotFoundError): return
    raise AssertionError('Altered workflow accepted')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True); p.add_argument('--receipt',type=Path,required=True)
    a = p.parse_args(); started=time.monotonic(); root=a.output.resolve(); root.mkdir(exist_ok=False)
    future_path=Path('metadata/baliphy_joint_node_logger_future_models_20261003_v4.json')
    future=json.loads(future_path.read_text()); verify(future)
    original_path=Path(future['source_plan']); original=json.loads(original_path.read_text()); verify(original)
    roles=json.loads(Path(future['future_roles']).read_text())
    originals=json.loads(Path(original['jobs']).read_text()); jobs=build_jobs(roles,originals)
    design_negatives=[]
    for case in ['missing_role','duplicate_id','duplicate_source','duplicate_seed','old_seed','bad_seed',
                 'changed_prior','changed_alias','changed_role','changed_program','invented_execution']:
        changed=copy.deepcopy(roles)
        if case=='missing_role': changed.pop()
        elif case=='duplicate_id': changed[1]['chain']['chain_id']=changed[0]['chain']['chain_id']
        elif case=='duplicate_source': changed[1]['source_chain_id']=changed[0]['source_chain_id']
        elif case=='duplicate_seed': changed[1]['chain']['seed']=changed[0]['chain']['seed']
        elif case=='old_seed': changed[0]['chain']['seed']=originals[0]['chain']['seed']
        elif case=='bad_seed': changed[0]['chain']['seed']=True
        elif case=='changed_prior': changed[0]['chain']['prior_label']='invented'
        elif case=='changed_alias': changed[0]['chain']['original_configuration_ids']=[]
        elif case=='changed_role': changed[0]['chain']['chain']=5
        elif case=='changed_program': changed[0]['chain']['program_sha256']='0'*64
        else: changed[0]['native_execution_launched']=True
        rejected(lambda:build_jobs(changed,originals)); design_negatives.append(case)
    gate_path=Path('metadata/baliphy_joint_node_logger_software_validation_20261003_v4.json')
    gate=json.loads(gate_path.read_text()); verify({'pins':gate['source_hashes']})
    fixture=Path('data/software_audits/baliphy-joint-node-logger-20261003-v4').resolve()
    mapping=root/'fixture_mapping.tsv'
    with mapping.open('w') as handle:
        writer=csv.DictWriter(handle,delimiter='\t',fieldnames=['guide','family','dataset','level','source_node','retained_set_json'])
        writer.writeheader()
        for i,n in enumerate([2,3,4,5]):
            writer.writerow(dict(guide='profile',family='software',dataset='whole',level=i,
                source_node='n'+str(i),retained_set_json=json.dumps(list('abcde')[:n])))
    bindings={str(path):sha(path) for path in [Path(__file__),future_path,original_path,Path(future['future_roles']),
        Path(original['jobs']),gate_path,mapping,fixture/'alignment.faa',fixture/'tree.nwk']}
    fixture_rows={}; fixture_arrays={}; native_negatives=[]
    for index,prior in enumerate(['broad','centered','package']):
        receipt=fixture/'native'/prior/'joint_node_logger'/'attempt-0001'/'receipt.json'
        config=json.loads((receipt.parent.parent/'configuration.json').read_text())
        program=fixture/(prior+'-joint_node_logger.hs')
        base=next(j['chain'] for j in jobs if j['chain']['prior_label']==prior)
        chain=dict(base,chain_id='software-'+prior,seed=20267001+index,family='software',proteins=5,
            original_configuration_ids=['software-whole-'+prior],effective_input_group='software',
            alignment=str(fixture/'alignment.faa'),alignment_sha256=sha(fixture/'alignment.faa'),
            tree=str(fixture/'tree.nwk'),tree_sha256=sha(fixture/'tree.nwk'),
            program=str(program),program_sha256=sha(program))
        job=dict(chain=chain,source_seed=20266001+index,memory_reservation_bytes=12*2**30,config=config)
        export=root/'actual_fixtures'/prior
        row=inspect(job,receipt,'software',mapping,export)
        assert row['status']==SUCCESS and len(row['joint_frames'])==3
        assert inspect(job,receipt,'software',mapping,export)==row
        fixture_rows[prior]=row; fixture_arrays[prior]=[]
        for frame in row['joint_frames']:
            with np.load(frame['projection_array'],allow_pickle=False) as arrays:
                fixture_arrays[prior].append({key:arrays[key].copy() for key in arrays.files})
        bindings[str(receipt)]=sha(receipt); bindings[str(receipt.parent.parent/'configuration.json')]=sha(receipt.parent.parent/'configuration.json')
        for name,h in json.loads(receipt.read_text())['artifacts'].items(): bindings[str(receipt.parent/name)]=h
        # Native artifacts remain unchanged: decoder errors are explicitly
        # mocked here, and never reported as actual native crashes.
        with patch('baliphy_joint_sampler_qualification.joint_frames',side_effect=ValueError('Synthetic malformed joint JSON')):
            bad=inspect(job,receipt,'software',mapping,root/'invalid_fixture'/prior)
            assert bad['status']==INVALID and bad['joint_frames']==[]
            assert bad['legacy_output_status']=='full_short_sampler_output_integrity_checked_not_posterior'
            assert not (root/'invalid_fixture'/prior).exists()
        native_negatives.append(prior+':invalid_joint_record_retained')
        altered=copy.deepcopy(row); altered['joint_frames'][0]['projection_array_sha256']='0'*64
        array_file=export/'frame-0.npz'; original_bytes=array_file.read_bytes()
        with array_file.open('wb') as handle: np.savez(handle,**{k:(v.astype(np.int16) if k=='categories' else v) for k,v in fixture_arrays[prior][0].items()})
        rejected(lambda:inspect(job,receipt,'software',mapping,export))
        array_file.write_bytes(original_bytes)
        native_negatives.append(prior+':changed_export_dtype_rejected')
    # Full real metadata through all producer/reader paths, with explicit mocks
    # for both native execution and output admission. These are not 1620 native
    # fixture runs or accepted fungal joint samples.
    job_path=root/'jobs.json'; write_json(job_path,jobs)
    output=root/'full_mock_grid'; plan_path=root/'mock_plan.json'
    plan=dict(jobs=str(job_path),output=str(output),mapping=str(mapping),pins={},scope='Explicit full-grid software mocks only',
        resources=dict(cpus=16,memory_gib=200,reservation_capacity_gib=192,workers=16,minimum_free_disk_gib=1))
    write_json(plan_path,plan); failures={jobs[0]['chain']['chain_id'],jobs[4]['chain']['chain_id']}
    def mock_inspect(job,receipt,digest,mapping,export):
        chain=job['chain']; prior=chain['prior_label']; row=copy.deepcopy(fixture_rows[prior])
        row.update(chain_id=chain['chain_id'],effective_input_group=chain['effective_input_group'],
            model_input_identity=chain['effective_input_group']+'-'+prior,chain_role=chain['chain'],
            family=chain['family'],original_configuration_ids=chain['original_configuration_ids'],
            seed=chain['seed'],source_seed=job['source_seed'],memory_reservation_bytes=job['memory_reservation_bytes'],plan_sha256=digest)
        if chain['chain_id'] in failures:
            row.update(status='unsuccessful_sampler_qualification_attempt_retained',exit_code=1,
                saved_alignments=0,candidate_frames=0,joint_frames=[],ancestral_categories_available=False,
                same_record_sequence_category_correspondence=False,legacy_output_status='unsuccessful_sampler_qualification_attempt_retained')
            return row
        export=Path(export); existing=export.exists()
        if not existing: export.mkdir(parents=True)
        for index,frame in enumerate(row['joint_frames']):
            path=export/('frame-'+str(frame['iteration'])+'.npz')
            if not existing: write_arrays(path,fixture_arrays[prior][index])
            verify_arrays(path,fixture_arrays[prior][index])
            frame.update(projection_array=str(path),projection_array_sha256=sha(path))
        return row
    def mock_attempt(folder,config):
        job=next(j for j in jobs if j['chain']['chain_id']==Path(folder).name)
        return Path(fixture_rows[job['chain']['prior_label']]['native_receipt'])
    serialization=[]
    with patch.object(workflow,'startup_gate',return_value={}),patch.object(workflow,'runtime_caps',return_value={}),\
            patch.object(workflow,'inspect',side_effect=mock_inspect),patch.object(workflow,'run_attempt',side_effect=mock_attempt):
        producer=workflow.run(plan_path); reader=workflow.run(plan_path,reader=True)
        assert producer['checked_sampler_attempts']==1618 and producer['unsuccessful_sampler_attempts']==2
        assert producer['complete_quartets']==403 and producer['unresolved_quartets']==2
        assert producer['joint_saved_frames']==4854
        rejected(lambda:workflow.run(plan_path)); serialization.append('producer_restart')
        rejected(lambda:workflow.run(plan_path,reader=True)); serialization.append('reader_restart')
        # Preserve actual successful reader; allow controlled re-entry solely
        # within this synthetic audit to check malformed role/array rejection.
        reader_file=output/'readback.json'; reader_bytes=reader_file.read_bytes(); reader_file.unlink()
        checkpoint=output/'chains'/(jobs[1]['chain']['chain_id']+'.json'); checkpoint_bytes=checkpoint.read_bytes()
        checkpoint.unlink(); rejected(lambda:workflow.run(plan_path,reader=True)); checkpoint.write_bytes(checkpoint_bytes)
        serialization.append('missing_checkpoint')
        changed=json.loads(checkpoint_bytes); changed['seed']+=1; checkpoint.write_text(json.dumps(changed))
        rejected(lambda:workflow.run(plan_path,reader=True)); checkpoint.write_bytes(checkpoint_bytes); serialization.append('changed_seed')
        changed=json.loads(checkpoint_bytes); changed['scientific_eligibility']=True; checkpoint.write_text(json.dumps(changed))
        rejected(lambda:workflow.run(plan_path,reader=True)); checkpoint.write_bytes(checkpoint_bytes); serialization.append('scientific_claim')
        reader_file.write_bytes(reader_bytes)
    # Closed-startup seed/metadata and full resource boundaries are exercised
    # with mocked closure admission, never fabricated completion journals.
    startup_rows=[dict(chain_id=j['chain']['chain_id'],source_seed=j['source_seed'],fresh_seed=j['chain']['seed'],
        chain_role=j['chain']['chain'],prior_label=j['chain']['prior_label'],effective_input_group=j['chain']['effective_input_group'],
        original_configuration_ids=j['chain']['original_configuration_ids'],scientific_eligibility=False,
        posterior_sampling_launched=False,status='reference_startup_homology_density_and_representation_checked') for j in jobs]
    startup_root=root/'mock_startup'; startup_root.mkdir(); write_json(startup_root/'dispositions.json',startup_rows)
    observer_root=root/'mock_observer'; observer_root.mkdir(); observations=observer_root/'observations.jsonl'
    snapshot=dict(sequence=0,cgroup=dict(limits={'cpu.max':'1600000 100000','memory.max':str(200*2**30),'memory.swap.max':'0'},
        memory_bytes={'memory.swap.current':0},memory_events={key:0 for key in ['max','oom','oom_kill','oom_group_kill']}))
    observations.write_text(json.dumps(snapshot)+'\n')
    startup_closed=dict(full_chains=1620,full_quartets=405,effective_inputs=135,original_configuration_aliases=324,
        validated_startups=1620,unsuccessful_startups=0,complete_startup_quartets=405,unresolved_startup_quartets=0,posterior_sampling_launched=False)
    sampler_closed=dict(full_chains=1620,full_quartets=405,effective_inputs=135,original_configuration_aliases=324,
        checked_sampler_attempts=1618,unsuccessful_sampler_attempts=2,posterior_qualified=False)
    observer_closed=dict(expected_roles=1620,roles_with_native_attempt=1620,roles_without_native_attempt=0,
        roles_with_live_observation=1619,roles_without_live_observation=1,posterior_qualified=False,
        sampler_controller_terminal_status='verified_original_terminal_success_with_bound_completed_artifacts',
        maximum_observed_cgroup_reported_peak_bytes=16*2**30)
    replay_closed=dict(full_chains=1620,checked_chains=1618,unresolved_chains=2,full_groups=405,effective_inputs=135,
        original_configuration_aliases=324,posterior_qualified=False,ancestral_categories_available=False,tip_logs_joint_trajectory_available=False)
    stage_plan=dict(jobs=str(job_path),startup_plan='startup',historical_sampler_plan='sampler',historical_resource_plan='observer',historical_replay_plan='replay')
    closures={'startup':(dict(output=str(startup_root)),startup_closed,{str(startup_root/'dispositions.json'):sha(startup_root/'dispositions.json')}),
        'sampler':({},sampler_closed,{}),
        'observer':(dict(output=str(observer_root),sampler_plan='sampler'),observer_closed,{str(observations):sha(observations)}),
        'replay':(dict(sampler_plan='sampler'),replay_closed,{})}
    gate_negatives=[]
    with patch('baliphy_joint_sampler_gates.closed_stage',side_effect=lambda path,status:closures[path]):
        prerequisites(stage_plan)
        for case in ['unsuccessful_startup','missing_sampler_role','missing_resource_attempt','replay_scientific_claim']:
            if case=='unsuccessful_startup': record,key,value=startup_closed,'unsuccessful_startups',1
            elif case=='missing_sampler_role': record,key,value=sampler_closed,'checked_sampler_attempts',1617
            elif case=='missing_resource_attempt': record,key,value=observer_closed,'roles_with_native_attempt',1619
            else: record,key,value=replay_closed,'posterior_qualified',True
            old=record[key]; record[key]=value; rejected(lambda:prerequisites(stage_plan)); record[key]=old; gate_negatives.append(case)
    for case in ['observed_oom','observed_swap','wrong_group_cap']:
        bad=copy.deepcopy(snapshot)
        if case=='observed_oom': bad['cgroup']['memory_events']['oom']=1
        elif case=='observed_swap': bad['cgroup']['memory_bytes']['memory.swap.current']=1
        else: bad['cgroup']['limits']['memory.max']='1'
        observations.write_text(json.dumps(bad)+'\n')
        rejected(lambda:resource_boundaries(dict(output=str(observer_root)),observer_closed,{str(observations):sha(observations)}))
        gate_negatives.append(case)
    observations.write_text(json.dumps(snapshot)+'\n')
    for path in ['scripts/baliphy_joint_sampler_qualification.py','scripts/baliphy_joint_sampler_gates.py',
        'scripts/run_baliphy_joint_sampler_qualification.py','scripts/reference_sampler_memory_budget.py',
        'scripts/independent_joint_ancestral_frames.py','scripts/baliphy_reference_sampler_qualification.py']:
        bindings[path]=sha(path)
    for path in [job_path,plan_path,output/'receipt.json',output/'readback.json']:
        bindings[str(path)]=sha(path)
    usage=resource.getrusage(resource.RUSAGE_SELF)
    result=dict(status='passed_full_joint_sampler_qualification_software_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        source_hashes=bindings,full_roles=1620,full_quartets=405,full_source_configurations_checked=True,
        actual_retained_native_fixture_roles=3,actual_retained_native_fixture_frames=9,
        full_mock_producer_reader_serialization_checked=True,artificial_failed_roles_retained=2,
        artificial_complete_quartets=403,artificial_unresolved_quartets=2,
        altered_designs_rejected=design_negatives,actual_fixture_decoder_export_negatives=native_negatives,
        serialization_cases_rejected=serialization,prerequisite_gate_negatives=gate_negatives,
        elapsed_seconds=time.monotonic()-started,self_cpu_seconds=usage.ru_utime+usage.ru_stime,
        self_peak_rss_bytes=usage.ru_maxrss*1024,new_native_inference_runs=0,posterior_qualified=False,scientific_eligibility=False,
        scope='All1620realmetadata/configs source inverse, disjointseed and resource classes;3retainedactualnative20iterationfixtures and9jointframes; full producer/readback with explicitly mocked native execution/admission and prerequisite closures;2artificialfailedroles retained. Frozen existing jobs unchanged. No full native joint production, actual new closure journals, posterior mixing/likelihood/model acceptance, GPU or charges.')
    with a.receipt.open('x') as handle: handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)


if __name__=='__main__': main()
