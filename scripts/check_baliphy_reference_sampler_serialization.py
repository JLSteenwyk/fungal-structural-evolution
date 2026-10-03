#!/usr/bin/env python3
"""Full serialized-reader and startup-gate contracts using labeled artificial rows."""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
from unittest.mock import patch

from ancestral_chain_attempt import sha, write_json
from baliphy_reference_sampler_qualification import summarize
import run_baliphy_reference_sampler_qualification as workflow
from run_baliphy_reference_preflight import verify


def rejected(action):
    try:action()
    except (AssertionError,KeyError,ValueError):return
    raise AssertionError('Altered qualification export or startup gate accepted')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--native-validation',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();native=json.loads(args.native_validation.read_text());verify(dict(pins=native['source_hashes']))
    assert native['status']=='passed_full_reference_sampler_qualification_software_contracts'
    root=args.output.resolve();root.mkdir(exist_ok=False)
    startup_path=Path('metadata/baliphy_reference_preflight_plan_20261003.json');startup=json.loads(startup_path.read_text())
    old_jobs=json.loads(Path(startup['jobs']).read_text())
    jobs=[dict(chain=dict(j['chain'],seed=j['fresh_seed']),source_seed=j['chain']['seed'],
        memory_reservation_bytes=(48 if j['chain']['family'] in ['OG0000972','OG0001082'] else 12)*2**30) for j in old_jobs]
    jobs_path=root/'jobs.json';write_json(jobs_path,jobs)
    stage=root/'stage';stage.mkdir();(stage/'chains').mkdir()
    plan_path=root/'plan.json';plan=dict(output=str(stage),jobs=str(jobs_path),mapping='Software mapping is mocked',pins={},
        resources=dict(workers=16,reservation_capacity_gib=192),scope='Artificial full-grid serialization fixture; native reconstruction mocked')
    write_json(plan_path,plan);digest=sha(plan_path);rows=[]
    for i,j in enumerate(jobs):
        c=j['chain'];template=next(r for r in native['native_rows'] if r['prior_label']==c['prior_label'])
        row=copy.deepcopy(template);row.update(chain_id=c['chain_id'],effective_input_group=c['effective_input_group'],
            model_input_identity=c['effective_input_group']+'-'+c['prior_label'],chain_role=c['chain'],family=c['family'],
            original_configuration_ids=c['original_configuration_ids'],seed=c['seed'],source_seed=j['source_seed'],
            memory_reservation_bytes=j['memory_reservation_bytes'],plan_sha256=digest)
        if i in [0,4]:
            row.update(status='unsuccessful_sampler_qualification_attempt_retained',exit_code=1,saved_alignments=0,candidate_frames=0)
            row.pop('sample_audit',None)
        rows.append(row)
    rows.sort(key=lambda r:r['chain_id']);by_id={r['chain_id']:r for r in rows};summary=summarize(rows)
    assert (summary['checked_sampler_attempts'],summary['unsuccessful_sampler_attempts'],summary['complete_quartets'])==(1618,2,403)
    for r in rows:write_json(stage/'chains'/(r['chain_id']+'.json'),r)
    fixture_root=Path(native['native_rows'][0]['native_receipt']).parents[4]
    events=json.loads((fixture_root/'full_scheduler_stress_events.json').read_text())
    ledger=workflow.verify_ledger(events,jobs,192*2**30,16)
    journal=stage/'memory_reservations.jsonl';journal.write_text(''.join(json.dumps(e)+'\n' for e in events))
    exported=stage/'dispositions.json';write_json(exported,rows);write_json(stage/'stage_plan.json',plan)
    producer=dict(status=workflow.PRODUCER_STATUS,plan_sha256=digest,**summary,reservation_audit=ledger,
        artifacts={str(p.relative_to(stage)):sha(p) for p in [exported,stage/'stage_plan.json',journal,*sorted((stage/'chains').glob('*.json'))]},
        scientific_eligibility=False)
    write_json(stage/'receipt.json',producer);negatives=[]
    def inspect_fixture(job,*ignored):return by_id[job['chain']['chain_id']]
    with patch.object(workflow,'runtime_caps',return_value={'software_fixture_only':True}),patch.object(
        workflow,'startup_gate',return_value={}),patch.object(workflow,'inspect',side_effect=inspect_fixture):
        original=exported.read_bytes()
        for name in ['missing_role','duplicated_role','changed_seed','removed_failed_role','changed_alias','scientific_acceptance','changed_memory_class']:
            bad=copy.deepcopy(rows)
            if name=='missing_role':bad.pop()
            elif name=='duplicated_role':bad[1]=copy.deepcopy(bad[0])
            elif name=='changed_seed':bad[1]['seed']+=1
            elif name=='removed_failed_role':bad[0]['status']='full_short_sampler_output_integrity_checked_not_posterior'
            elif name=='changed_alias':bad[1]['original_configuration_ids']=['invented']
            elif name=='scientific_acceptance':bad[1]['scientific_eligibility']=True
            else:bad[1]['memory_reservation_bytes']=48*2**30
            write_json(exported,bad);rejected(lambda:workflow.run(plan_path,reader=True));negatives.append(name)
            exported.write_bytes(original)
        reader=workflow.run(plan_path,reader=True);assert all(reader[k]==v for k,v in summary.items())
        rejected(lambda:workflow.run(plan_path,reader=True));negatives.append('completed_reader_restart')
    # Test the real startup-gate implementation separately. Its files and
    # journal entries are artificial; actual prerequisite closure remains
    # a production requirement, never inferred from this software fixture.
    gate_root=root/'startup';gate_root.mkdir();dispositions=gate_root/'dispositions.json'
    startup_rows=[dict(chain_id=j['chain']['chain_id'],fresh_seed=j['chain']['seed'],source_seed=j['source_seed'],
        effective_input_group=j['chain']['effective_input_group'],prior_label=j['chain']['prior_label'],
        chain_role=j['chain']['chain'],original_configuration_ids=j['chain']['original_configuration_ids'],
        status='reference_startup_homology_density_and_representation_checked') for j in jobs]
    write_json(dispositions,startup_rows)
    archive=root/'startup_archive.json';proof=dict(services=[{'software_fixture':True},{'software_fixture':True}],source_hashes={str(dispositions):sha(dispositions)})
    write_json(archive,proof)
    completion=root/'startup_completed.json';closed=dict(status='complete_verified_full_reference_startup_footer_replay',
        full_chains=1620,full_quartets=405,effective_inputs=135,original_configuration_aliases=324,validated_startups=1620,
        unsuccessful_startups=0,complete_startup_quartets=405,unresolved_startup_quartets=0,
        full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive))
    write_json(completion,closed);footer_plan=root/'startup_plan.json'
    write_json(footer_plan,dict(output=str(gate_root),completion=str(completion),pins={}))
    gate_plan=dict(plan,startup_plan=str(footer_plan));workflow.startup_gate(gate_plan)
    for name in ['incomplete_startup_count','failed_startup_count','incomplete_quartet','wrong_closure_status']:
        bad=dict(closed)
        if name=='incomplete_startup_count':bad['validated_startups']=1619
        elif name=='failed_startup_count':bad['unsuccessful_startups']=1
        elif name=='incomplete_quartet':bad['complete_startup_quartets']=404
        else:bad['status']='invented_completion'
        write_json(completion,bad);rejected(lambda:workflow.startup_gate(gate_plan));negatives.append(name)
    write_json(completion,closed)
    for name in ['changed_startup_seed','invalid_startup_disposition']:
        bad=copy.deepcopy(startup_rows)
        if name=='changed_startup_seed':bad[0]['fresh_seed']+=1
        else:bad[0]['status']='invalid_native_startup_retained'
        write_json(dispositions,bad);proof['source_hashes'][str(dispositions)]=sha(dispositions);write_json(archive,proof)
        closed['full_hash_archive_sha256']=sha(archive);write_json(completion,closed)
        rejected(lambda:workflow.startup_gate(gate_plan));negatives.append(name)
    write_json(dispositions,startup_rows);proof['source_hashes'][str(dispositions)]=sha(dispositions);write_json(archive,proof)
    closed['full_hash_archive_sha256']=sha(archive);write_json(completion,closed);workflow.startup_gate(gate_plan)
    combined=dict(native);combined.update(checked_utc=datetime.now(timezone.utc).isoformat(),
        full_serialized_reader_checked=True,full_startup_admission_gate_checked=True,
        serialization_and_prerequisite_alterations_rejected=negatives,
        native_validation=str(args.native_validation),native_validation_sha256=sha(args.native_validation),
        serialization_scope='Full1620real metadata identities with explicitly artificial status rows and mock reconstruction. Real reader export/hash/reservation accounting and startup gate implementation exercised, including two retained failures and14 rejection cases. This is not actual full-grid native output or journal closure.')
    combined['source_hashes']=dict(native['source_hashes'])
    for p in [Path(__file__),args.native_validation,startup_path,*sorted(root.rglob('*.json')),journal]:combined['source_hashes'][str(p)]=sha(p)
    with args.receipt.open('x') as f:f.write(json.dumps(combined,indent=2)+'\n')
    print(json.dumps(dict(status=combined['status'],full_serialized_reader_checked=True,
        full_startup_admission_gate_checked=True,rejection_cases=negatives)),flush=True)


if __name__=='__main__':main()
