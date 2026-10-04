#!/usr/bin/env python3
"""Freeze full-scope parallel numerical plans after original software/source closure."""
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import psutil

from ancestral_chain_attempt import sha
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind, verify


def write(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2);f.write('\n')


def main():
    oldpath=Path('metadata/full_weighted_covariance_qualification_plan_20261004_v1.json')
    old=json.loads(oldpath.read_text());verify(old['pins'])
    software=[];pins=dict(old['pins']);bind(pins,oldpath)
    for stem,version,session in [('weighted_parallel_qualification',2,48622),
            ('weighted_parallel_source_guards',1,65804)]:
        paths=[Path('metadata/'+stem+'_software_'+kind+'_20261004_v'+str(version)+'.json') for kind in ['validation','execution','transport']]
        v,e,t=[json.loads(p.read_text()) for p in paths]
        assert e['status']=='exited_zero_with_receipt' and e['exit_code']==0 and not e['timed_out']
        assert sha(paths[0])==e['receipt_sha256'] and t['original_tool_terminal_exit_code']==0
        assert t['original_tool_session_id']==session and t['wrapper']==e['wrapper']
        assert t['invocation_id']==e['invocation_id'] and t['whole_wrapper_initial_and_terminal_payloads_matched']
        assert t['manager_start_records']==t['manager_completion_records']==1
        for r in [v,e,t]:verify(r['source_hashes']);verify(r.get('artifacts',{}))
        for path in paths:bind(pins,path)
        software.append(v)
    q,g=software
    assert q['status']=='passed_complete_parallel_four_control_software_qualification_v1'
    assert (q['serial_parallel_identical_audit_records'],q['serial_parallel_identical_setting_links'])==(18000,72000)
    assert q['independent_latent_readbacks']==3 and q['qualified_basis_dimensions']==[4,5,6]
    assert len(q['rehashed_private_corruptions_rejected'])==19 and q['retained_original_numeric_failure_captures']==13
    assert q['completed_stage_restart_refusals']==6 and q['unchanged_comparison_rtol']==3e-9 and q['unchanged_comparison_atol']==2e-8
    assert g['status']=='passed_private_fork_cached_source_mutation_guards_v1'
    assert [(c['case'],c['rejected']) for c in g['cases']]==[('control',False),('before',True),('after',True)]
    assert g['parent_numeric_source_preserved'] and g['original_input_files_preserved']
    sourcepath=Path(old['source_census_completion']);source=json.loads(sourcepath.read_text())
    assert source['status']=='complete_verified_full_four_control_covariance_source_census_v1'
    assert source['cohorts']==4340 and source['exact_process_journals_checked']==2
    archivepath=Path(source['full_hash_archive']);assert sha(archivepath)==source['full_hash_archive_sha256']
    archive=json.loads(archivepath.read_text());verify(archive['source_hashes'])
    assert len(archive['source_hashes'])==source['bound_source_hashes']==13119
    bind(pins,sourcepath);bind(pins,archivepath)
    modules=project_sources(pins,[Path(__file__),
        Path('scripts/full_weighted_covariance_qualification_parallel_v1.py'),
        Path('scripts/close_full_weighted_parallel_stage_v1.py'),
        Path('scripts/record_full_weighted_parallel_launch_v1.py')])
    available=psutil.virtual_memory().available;disk=psutil.disk_usage('.').free
    assert available>=200*2**30 and disk>=256*2**30
    resources=dict(checked_utc=datetime.now(timezone.utc).isoformat(),cpus=16,memory_gib=200,
        workers=16,worker_address_space_gib=12,reservation_capacity_gib=192,parent_headroom_gib=8,
        parent_soft_address_space_gib=8,parent_hard_address_space_gib=24,swap_gib=0,blas_threads=1,
        cpu_seconds_per_process=604800,per_file_limit_mib=512,minimum_free_disk_gib=256,
        minimum_available_ram_gib=200,output_and_scratch_allowance_gib=160,
        available_ram_gib=available/2**30,free_disk_gib=disk/2**30,
        uncompressed_audit_budget_bytes=5208000*16384,uncompressed_link_budget_bytes=24883200*2048,
        original_serial_process_vm_gib=psutil.Process(3078625).memory_info().vms/2**30,
        original_serial_process_rss_gib=psutil.Process(3078625).memory_info().rss/2**30,
        estimated_core_hours_per_stage=[30,100],planning_wall_hours_per_stage=[3,12],
        full_parallel_scaling_measured=False,gpu=False,new_cost_usd=0,
        scope='Sixteen fork workers, each limited to12GiB address space, plus8GiB parent '
              'soft address space and200GiB cgroup/no swap. Parent hard limit24GiB permits '
              'children to install12GiB limits; parent remains limited to8GiB. BLAS1. '
              'Full original scope and output allowance unchanged. Runtime range is a '
              'planning estimate informed by original serial progress, not a calibrated ETA. '
              'The old original serial job and its queued timing remain unchanged.')
    rp=Path('metadata/full_weighted_parallel_qualification_resources_20261004_v1.json');write(rp,resources);bind(pins,rp)
    output='results/phylogeny/full-four-control-covariance-parallel-qualification-20261004-v1'
    assert not Path(output).exists()
    pp=Path('metadata/full_weighted_parallel_qualification_plan_20261004_v1.json')
    scope=old['scope']+' Parallel scheduling only: every original cohort, audit, setting and guard retained. '
    scope+='Cached numeric source hashes checked before/after every cohort; first failure stops new admissions, no automatic retry. '
    scope+='Independent reader rebuilds all numerical cases and binds all producer/reader checkpoints before original two-journal closure. '
    scope+='Original failed timing/diagnostics and active serial numerical/timing stages are preserved. All eight aims remain incomplete.'
    plan=dict(old,output=output,resources=resources,pins=pins,scope=scope,
        completion='metadata/full_weighted_parallel_qualification_completed_20261004_v1.json',
        closure_script='scripts/close_full_weighted_parallel_stage_v1.py',
        software_qualification='metadata/weighted_parallel_qualification_software_transport_20261004_v2.json',
        source_guard_qualification='metadata/weighted_parallel_source_guards_software_transport_20261004_v1.json',
        source_modules=sorted(map(str,modules)))
    verify(pins);write(pp,plan)
    inventory=json.loads(Path('metadata/full_weighted_covariance_source_census_launches_20261004_v1.json').read_text())
    dependency=inventory['launches'][2]
    command=['/usr/bin/prlimit','--as='+str(8*2**30)+':'+str(24*2**30),
        '--cpu=604800','--fsize='+str(512*2**20),'--',sys.executable,
        'scripts/full_weighted_covariance_qualification_parallel_v1.py','--plan',str(pp)]
    wp=Path('metadata/full_weighted_parallel_producer_wait_plan_20261004_v1.json')
    write(wp,dict(dependencies=[dependency],command=command,
        pins={str(p):sha(p) for p in [pp,Path(dependency),Path('scripts/run_after_verified_dependencies_v2.py')]}))
    print(json.dumps(dict(status='prepared_full_parallel_numerical_producer_not_yet_launched',
        plan=str(pp),plan_sha256=sha(pp),wait_plan=str(wp),source_bindings_rechecked=13119,
        expected=plan['expected'],resources=resources),indent=2))


if __name__=='__main__':main()
