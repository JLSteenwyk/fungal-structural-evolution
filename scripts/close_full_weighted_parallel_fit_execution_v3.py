#!/usr/bin/env python3
"""Close complete fit/readback output and original bounded role custody together."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from ancestral_chain_attempt import sha
from full_weighted_fit_exports import atomic, PRODUCER, READER, SUMMARY
from reference_measurement_union_sources import bind, verify
from weighted_parallel_fit_execution_runtime_v3 import configuration, native_limits, role_command, validate_native_receipt


def native_custody(folder,role,command,limits,resources,inputs,bindings):
    """Bind an original guard attempt; also usable for non-fit boundary checks."""
    rp=Path(folder)/'receipt.json';record=json.loads(rp.read_text())
    assert record['status']=='original_native_weighted_fit_role_exited_zero' and record['exit_code']==0
    assert record['role']==role and record['scientific_eligibility'] is record['restarted'] is False
    assert record['child']['created'] is not None and record['child']['pid']>0
    config=json.loads((rp.parent/'configuration.json').read_text());process=json.loads((rp.parent/'process.json').read_text())
    validate_native_receipt(record,config,process,command,limits,resources)
    inputs=dict(inputs)
    for p in [sys.executable,'/usr/bin/prlimit']:bind(inputs,p)
    assert record['source_hashes']==inputs
    artifacts=[rp.parent/name for name in ['configuration.json','process.json','stdout.log','stderr.log']]
    assert record['artifacts']=={str(p):sha(p) for p in artifacts}
    for p,d in record['source_hashes'].items():bind(bindings,p,d)
    for p,d in record['artifacts'].items():bind(bindings,p,d)
    bind(bindings,rp)
    verify(bindings)
    return record


def role_receipt(folder,role,operation,fit,resources,bindings):
    expected=native_limits(resources,role)
    command=['/usr/bin/prlimit','--as='+str(expected['address_space_bytes']),'--cpu='+str(expected['cpu_seconds']),
        '--fsize='+str(expected['per_file_bytes']),'--',*role_command(operation['fit_plan'],fit,role)]
    inputs=dict(bindings)
    if role=='reader':bind(inputs,Path(fit['output'])/'receipt.json')
    directory=Path(folder)/role
    record=native_custody(directory,role,command,expected,resources,inputs,bindings)
    output_path=directory/'role_output.json';output=json.loads(output_path.read_text())
    target=Path(fit['output'])/('receipt.json' if role=='producer' else 'readback.json')
    assert output==dict(receipt=str(target),sha256=sha(target),status=PRODUCER if role=='producer' else READER,
        fit_contract=operation['fit_contract'],scientific_eligibility=False)
    bind(bindings,output_path);verify(bindings)
    return record


def run(plan_path):
    plan,operation,fit,resources,bindings=configuration(plan_path)
    assert not Path(plan['completion']).exists()
    original_bindings=dict(bindings)
    for role in ['producer','reader']:
        local=dict(original_bindings)
        record=role_receipt(plan['runtime_root'],role,operation,fit,resources,local)
        assert record['wrapper']['cmdline']==[sys.executable,'scripts/weighted_parallel_fit_execution_runtime_v3.py','--plan',str(plan_path),'--role',role]
        for p,d in local.items():bind(bindings,p,d)
    inventory=json.loads(Path(plan['role_launch_inventory']).read_text())
    assert inventory['controller_plan_sha256']==sha(plan_path) and len(inventory['launches'])==2
    # The dependency wrappers' original journals are checked by the generic
    # closer; native identities/caps/logs are additionally bound above.
    bind(bindings,plan['role_launch_inventory'])
    for launch in inventory['launches'][:2]:bind(bindings,launch)
    root=Path(fit['output']);producer=json.loads((root/'receipt.json').read_text());reader=json.loads((root/'readback.json').read_text())
    assert producer['source_contract']==reader['source_contract']==operation['fit_contract']
    assert reader['producer_receipt_sha256']==sha(root/'receipt.json')
    assert all(producer[k]==reader[k] for k in SUMMARY)
    for value in [producer,reader]:
        assert value['cached_numeric_inputs_preserved_for_all_cohorts'] is True and value['native_call_numeric_inputs_checked'] is True
        assert value['scientific_eligibility'] is value['nonuniform_weighting_accepted'] is value['component_variance_attribution_accepted'] is False
    marker=root/'independent_readback_completed.json'
    assert json.loads(marker.read_text())==dict(output=str(root/'readback.json'),output_sha256=sha(root/'readback.json'),
        plan_sha256=operation['fit_plan_sha256'],producer_receipt_sha256=sha(root/'receipt.json'))
    bind(bindings,marker)
    spec=dict(source_plan=operation['fit_plan'],producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status=PRODUCER,reader_status=READER,completed_status='complete_verified_full_four_control_shared_entity_fits_v1',
        summary_fields=SUMMARY,launches=inventory['launches'][:2],pins=bindings,output=plan['completion'],scope=plan['scope'])
    target=Path(plan['runtime_root'])/'full_fit_closure_plan.json';atomic(target,spec)
    subprocess.run([sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',str(target)],check=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
