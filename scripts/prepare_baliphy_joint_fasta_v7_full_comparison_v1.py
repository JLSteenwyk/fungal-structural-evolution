#!/usr/bin/env python3
"""Freeze a full failed-input candidate comparison, retaining all original native caps."""
from datetime import datetime, timezone
import json
from pathlib import Path

import psutil

from ancestral_chain_attempt import sha
from baliphy_joint_fasta_lines_v7 import reverse, transform
from reference_measurement_union_sources import bind, verify


def save(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def main():
    pp=Path('metadata/baliphy_joint_fasta_v7_full_comparison_plan_20261004_v1.json')
    rp=Path('metadata/baliphy_joint_fasta_v7_full_comparison_resources_20261004_v1.json')
    assert not pp.exists() and not rp.exists()
    gp=Path('metadata/baliphy_joint_fasta_v7_software_validation_20261004_v1.json')
    tp=Path('metadata/baliphy_joint_fasta_v7_software_transport_20261004_v1.json')
    gate,transport=[json.loads(p.read_text()) for p in [gp,tp]]
    assert gate['status']=='passed_reversible_rowwise_joint_fasta_v7_paired_native_controls'
    assert gate['full_programs_reversibly_checked']==405 and gate['native_prior_controls']==3
    assert transport['validation_sha256']==sha(gp) and transport['original_tool_terminal_exit_code']==0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    pins=dict(transport['source_hashes']);verify(pins)
    bp=Path('metadata/scalar_native_signal_debugger_plan_20261004_v1.json')
    baseline=json.loads(bp.read_text());verify(baseline['pins']);pins.update(baseline['pins'])
    original=Path(baseline['native_command'][baseline['native_command'].index('run')+1])
    prepared=Path('results/ancestral/prepared-joint-fasta-v7-full-comparison-20261004-v1')
    prepared.mkdir(exist_ok=False)
    model=(prepared/original.name).resolve();source=original.read_text();model.write_text(transform(source))
    assert reverse(model.read_text())==source
    command=list(baseline['native_command']);command[command.index('run')+1]=str(model)
    caps=dict(baseline['native_caps'],stack_soft_bytes=8*2**20)
    command=['/usr/bin/prlimit','--as='+str(caps['address_space_bytes']),
             '--cpu='+str(caps['cpu_seconds']),'--fsize='+str(caps['per_file_bytes']),
             '--stack='+str(caps['stack_soft_bytes'])+':unlimited','--',*command]
    root='results/ancestral/baliphy-joint-fasta-v7-full-input-comparison-20261004-v1'
    assert not Path(root).exists()
    scope=('One separate full622-tip failed-family comparison of row-wise joint FASTA construction. '
           'Original model, prior, inputs,seed1322110685,20iterations and48GiB AS/5674CPU-seconds/'
           '2GiB file/8MiB soft-unlimited hard stack retained; only pure logger construction changes. '
           'No debugger or phase marker is used in this comparison. No original attempt is retried '
           'or replaced. Successful20iterations alone are not a full24-failure repair, '
           'full-family sensitivity qualification, adequate posterior or biological result.')
    resources=dict(prepared_utc=datetime.now(timezone.utc).isoformat(),cpus=2,memory_gib=64,
        swap_gib=0,blas_threads=1,address_space_gib=64,cpu_seconds_per_stage=7200,
        wall_seconds_per_stage=7500,per_file_limit_mib=2048,
        minimum_available_ram_gib=64,minimum_free_disk_gib=256,
        available_ram_gib=psutil.virtual_memory().available/2**30,
        free_disk_gib=psutil.disk_usage('.').free/2**30,
        output_allowance_gib=4,estimated_wall_minutes=[2,15],estimate_uncalibrated=True,
        original_failed_role_wall_seconds=baseline['resources']['observed_original_role_wall_seconds'],
        maximum_native_cpu_seconds=caps['cpu_seconds'],new_cost_usd=0,gpu=False,scope=scope)
    assert resources['available_ram_gib']>=64 and resources['free_disk_gib']>=256
    save(rp,resources)
    for p in [Path(__file__),Path('scripts/baliphy_joint_fasta_lines_v7.py'),
              Path('scripts/run_baliphy_joint_fasta_v7_full_comparison_v1.py'),
              Path('scripts/ancestral_chain_attempt.py'),Path('scripts/read_baliphy_scalar_json_v6b.py'),
              Path('scripts/baliphy_joint_node_logger_v3.py'),Path('/usr/bin/prlimit'),
              original,model,gp,tp,bp,rp]:bind(pins,p)
    verify(pins)
    config=dict(command=command,timeout_seconds=7200,pins={str(Path(p).resolve()):h for p,h in pins.items()})
    plan=dict(status='prepared_candidate_full_input_comparison_after_native_controls',
        prepared_utc=datetime.now(timezone.utc).isoformat(),output=root,native_config=config,
        native_caps=caps,original_model=str(original),original_model_sha256=sha(original),
        diagnostic_model=str(model),diagnostic_model_sha256=sha(model),
        original_debugger_directory='results/ancestral/scalar-native-signal-debugger-20261004-v1/independent-chain-1',
        exact_model_reversibility_verified=True,resources=resources,pins=pins,scope=scope)
    save(pp,plan)
    print(json.dumps({k:v for k,v in plan.items() if k not in ['pins','resources','native_config']},indent=2))


if __name__=='__main__':main()
