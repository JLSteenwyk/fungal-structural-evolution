#!/usr/bin/env python3
"""Run the full CPU-only resource survey and queue complete readback/closure."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import psutil
from launch_full_triad_sequence_geometry_followup import create
from record_project_runtime_checkpoint_v4 import fingerprint,journal_terminal
from run_ortholog_pair_guide_comparison import sha


SUMMARY_FIELDS=['full_chains','effective_input_groups','full_quartets','selected_checked_chains','selected_failed_chains',
    'original_and_selected_recovery_attempts','complete_quartets','unresolved_quartets','proposed_fresh_chain_seeds',
    'proposed_iterations','observed_successful_elapsed_worker_seconds','observed_successful_native_alignment_bytes_size_only',
    'linear_successful_worker_seconds','linear_successful_native_alignment_bytes_size_only',
    'memory_requirements_for_two_failed_chains_unresolved','runtime_projection_uncalibrated','finish_eta','convergence_eta',
    'sampler_review_quartets','selected_successful_allocation_warning_chains','selected_successful_parameter_review_chains',
    'attempt_diagnostics','production_launch_allowed']


def launch(name,command,dependencies,source):
    prefix='metadata/'+name.replace('-','_');unit='fungal-'+name+'-20261003-v1.service'
    wait_path=prefix+'_wait_plan_20261003.json';launch_path=prefix+'_launch_20261003.json'
    assert not Path(launch_path).exists()
    create(wait_path,dict(dependencies=dependencies,command=command,
        pins={p:sha(p) for p in [source,*dependencies,'scripts/run_after_verified_dependencies_v2.py']}))
    cmd=['systemd-run','--user','--collect','--unit='+unit,'--working-directory='+str(Path.cwd()),
        '--property=CPUQuota=200%','--property=MemoryMax=8G','--property=MemorySwapMax=0',
        '--setenv=PYTHONUNBUFFERED=1','--setenv=OPENBLAS_NUM_THREADS=1','--setenv=OMP_NUM_THREADS=1','--setenv=MKL_NUM_THREADS=1',
        sys.executable,'scripts/run_after_verified_dependencies_v2.py','--plan',wait_path]
    subprocess.run(cmd,check=True);pid=0
    for _ in range(20):
        pid=int(subprocess.check_output(['systemctl','--user','show',unit,'-p','MainPID','--value'],text=True))
        if pid:break
        time.sleep(.1)
    assert pid>0;proc=psutil.Process(pid)
    expected=[sys.executable,'scripts/run_after_verified_dependencies_v2.py','--plan',wait_path];assert proc.cmdline()==expected
    cg=subprocess.check_output(['systemctl','--user','show',unit,'-p','ControlGroup','--value'],text=True).strip()
    limits={k:(Path('/sys/fs/cgroup')/cg.lstrip('/')/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
    assert limits=={'cpu.max':'200000 100000','memory.max':str(8*2**30),'memory.swap.max':'0'}
    create(launch_path,dict(unit=unit,pid=pid,created=proc.create_time(),cmdline=expected,plan=wait_path,
        plan_sha256=sha(wait_path),checked_utc=datetime.now(timezone.utc).isoformat(),actual_cgroup_limits=limits,systemd_launch_command=cmd))
    print(json.dumps(dict(launch=launch_path,pid=pid,unit=unit)),flush=True);return launch_path


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text())
    for path,d in plan['pins'].items():assert sha(path)==d,path
    grid=json.loads(Path(plan['grid_validation']).read_text())
    assert grid['status']=='passed_full_baliphy_horizon_resource_inventory_software_contracts'
    assert (grid['full_chains'],grid['full_quartets'],grid['attempts'],grid['fresh_seeds'])==(1620,405,1623,1620)
    assert grid['full_producer_and_serialized_reader_exercised'] and grid['completed_restarts_refused']
    assert len(grid['false_exports_and_malformed_traces_rejected'])==12
    for path,d in grid['source_hashes'].items():assert sha(path)==d,path
    assert plan['resources']['memory_gib']==8 and plan['resources']['cpu_per_stage']==2 and plan['proposed_iterations']==10000
    assert psutil.virtual_memory().available>=8*2**30 and shutil.disk_usage('.').free>=64*2**30
    assert not Path(plan['output']).exists()
    for path in plan['dependencies']:
        row=json.loads(Path(path).read_text());row['launch']=path;assert fingerprint(row) is None;journal_terminal(row)
    root=Path(plan['output'])
    producer=launch('baliphy-horizon-resources',[sys.executable,'scripts/prepare_baliphy_horizon_resources.py','--plan',str(a.plan)],plan['dependencies'],str(a.plan))
    reader=launch('baliphy-horizon-resources-readback',[sys.executable,'scripts/readback_baliphy_horizon_resources.py','--plan',str(a.plan)],[producer],str(a.plan))
    completion=dict(source_plan=str(a.plan),producer_receipt=str(root/'receipt.json'),independent_readback=str(root/'readback.json'),
        producer_status='complete_full_baliphy_horizon_resource_inventory_pending_readback',
        reader_status='passed_full_baliphy_horizon_resource_inventory_serialized_readback',
        completed_status='complete_verified_full_baliphy_horizon_resource_inventory',summary_fields=SUMMARY_FIELDS,
        launches=[producer,reader],pins={p:sha(p) for p in [str(a.plan),producer,reader,'scripts/close_full_triad_sequence_stage.py',
        'scripts/record_completed_process_handoffs_v2.py']},output=plan['completion'],scope=plan['scope'])
    create(plan['completion_plan'],completion)
    closer=launch('baliphy-horizon-resources-closure',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',plan['completion_plan']],[producer,reader],plan['completion_plan'])
    create(plan['launch_inventory'],dict(status='launched_full_baliphy_horizon_resource_inventory',
        source_plan=str(a.plan),source_plan_sha256=sha(a.plan),launches=[producer,reader,closer],
        native_sampling_launched=False,gpu=False,new_cost_usd=0))


if __name__=='__main__':main()
