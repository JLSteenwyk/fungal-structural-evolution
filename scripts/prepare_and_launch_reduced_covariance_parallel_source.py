#!/usr/bin/env python3
"""Queue the complete retained-basis grid after the alternate full source closure."""
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys

import psutil

from ancestral_chain_attempt import sha
from full_reduced_covariance_sources_v2 import SUMMARY_FIELDS
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from reference_measurement_union_sources import verify
from reduced_covariance_basis import validate_certificate
from full_covariance_qualification_sources import jsonl
from run_full_reduced_covariance_qualification_v2 import PRODUCER_STATUS,READER_STATUS


def main():
    gp=Path('metadata/reduced_covariance_parallel_source_software_validation_20261003_v3.json');gate=json.loads(gp.read_text())
    ep=Path('metadata/reduced_covariance_parallel_source_software_execution_20261003_v3.json');execution=json.loads(ep.read_text())
    assert gate['status']=='passed_full_exact_retained_covariance_parallel_source_v3_software_contracts'
    assert (gate['synthetic_pipeline_audits'],gate['synthetic_pipeline_setting_links'])==(600,1200)
    assert len(gate['alternate_source_cases_rejected'])==16
    for k in ['byte_exact_post_tampering_fixture_restoration','complete_positive_source_bindings_rehashed',
        'original_reader_path_not_required_or_written','inherited_numerical_envelopes_unchanged']:
        assert gate[k] is True
    verify(gate['source_hashes']);verify(execution['source_hashes']);verify(execution['artifacts'])
    assert execution['status']=='exited_zero_with_receipt' and execution['exit_code']==0 and execution['receipt_sha256']==sha(gp)
    # Preserve both actual waits. V2 passed at execution, but its post-mutation
    # fixture restoration finding prevents use as the admission gate.
    transports=[];bindings={}
    for v,session,inv in [(2,22774,'a7559ab6db91472eb95fbf28bce30296'),(3,35338,'89d26cba38224452b358f28b021bb7ff')]:
        xp=Path(f'metadata/reduced_covariance_parallel_source_software_execution_20261003_v{v}.json')
        x=json.loads(xp.read_text());unit=f'fungal-retained-covariance-parallel-source-software-20261003-v{v}.service'
        raw=subprocess.check_output(['journalctl','--user','-u',unit,'-o','json','--no-pager'],text=True)
        rows=[json.loads(l) for l in raw.splitlines()]
        exact=[r for r in rows if r.get('_PID')==str(x['wrapper']['pid']) and r.get('_CMDLINE')==' '.join(x['wrapper']['cmdline'])]
        assert exact and {r['_SYSTEMD_INVOCATION_ID'] for r in exact}=={inv}
        resources=[r for r in rows if r.get('USER_INVOCATION_ID')==inv and r.get('CPU_USAGE_NSEC')]
        assert resources and not any('Failed with result' in (r.get('MESSAGE') or '') or 'Main process exited' in (r.get('MESSAGE') or '') for r in rows if r.get('USER_INVOCATION_ID')==inv)
        jp=xp.with_suffix('')/'original-invocation-journal.jsonl';assert not jp.exists();jp.write_text(raw)
        transports.append(dict(version=v,actual_tool_terminal_exit_code=0,tool_session_id=session,invocation_id=inv,
            original_process_messages=len(exact),original_completion_resource_records=len(resources),
            execution=str(xp),execution_sha256=sha(xp),post_mutation_fixture_gate_accepted=v==3))
        bindings.update({str(xp):sha(xp),str(jp):sha(jp)})
    tp='metadata/reduced_covariance_parallel_source_software_transport_20261003_v3.json'
    create(tp,dict(status='verified_original_alternate_covariance_software_waits_with_preserved_restoration_finding',
        transports=transports,source_hashes=bindings,
        scope='Two actual original software tool waits and PID/create/command/invocation journals. V3 complete positive source graph rehashed; V2 test restoration finding retained. Borrowed execution wrapper has an inherited scope paragraph mentioning8680real certificates; these two checker commands themselves cover600synthetic audits/1200links. Full real certificate replay is separately performed by this preparation command, not attributed to those children. No production source edits or biological acceptance.'))
    primitive_path=Path('metadata/reduced_covariance_qualification_software_validation_20261003_v1.json')
    primitive=json.loads(primitive_path.read_text());verify(primitive['source_hashes']);verify(primitive['artifacts'])
    assert primitive['dense_nonnegative_cone_roundtrips']==80 and len(primitive['malformed_cases_rejected'])==16
    exact_completion='metadata/full_exact_covariance_folds_completed_20261003_v2.json';ec=json.loads(Path(exact_completion).read_text())
    assert (ec['cohorts'],ec['certificates'],ec['logical_cases'])==(4340,8680,75188)
    assert sha(ec['full_hash_archive'])==ec['full_hash_archive_sha256']
    exact_plan='metadata/full_exact_covariance_folds_plan_20261003_v2.json';exact=json.loads(Path(exact_plan).read_text())
    cert_path=Path(exact['output'])/'cohort_certificates.jsonl';archive=json.loads(Path(ec['full_hash_archive']).read_text())
    assert sha(cert_path)==archive['source_hashes'][str(cert_path)]
    counts={};total=0
    for c in jsonl(cert_path):
        n=len(validate_certificate(c));counts[str(n)]=counts.get(str(n),0)+1;total+=1
    assert total==8680 and counts=={'4':7232,'5':1448}
    dependency='metadata/parallel_covariance_readback_v1_closure_launch_20261003.json'
    resources=dict(checked_utc=datetime.now(timezone.utc).isoformat(),cpus=2,memory_gib=16,swap_gib=0,workers=1,
        address_space_gib=12,cpu_seconds=21600,per_file_limit_gib=1,blas_threads=1,output_allowance_gib=8,
        minimum_free_disk_gib=64,runtime_uncalibrated=True,finish_eta=None,new_cost_usd=0,gpu=False,
        available_memory_gib=psutil.virtual_memory().available/2**30,available_disk_gib=shutil.disk_usage('.').free/2**30,
        caveat='Full unchanged principal-Gram/envelope qualification after separate closed full parallel source proof. All1302000audits/6220800links retained. Source/proof dictionaries may dominate RAM; output8GiB is planning rather than a quota. No optimizer or accepted inference; existing V1 queue remains unchanged.')
    assert psutil.virtual_memory().available>=16*2**30 and shutil.disk_usage('.').free>=64*2**30
    own=['full_reduced_covariance_sources_v2','run_full_reduced_covariance_qualification_v2',
        'check_full_reduced_covariance_qualification_v3','prepare_and_launch_reduced_covariance_parallel_source',
        'reduced_covariance_basis','full_exact_covariance_sources','covariance_exact_folds_v2','covariance_basis_audit',
        'run_after_verified_dependencies_v2','close_full_triad_sequence_stage','record_completed_process_handoffs_v2']
    paths=[Path('scripts/'+n+'.py') for n in own]+[gp,ep,Path(tp),primitive_path,Path(exact_completion),Path(exact_plan),cert_path,
        Path('metadata/reduced_covariance_parallel_source_software_resources_20261003_v3.json'),
        Path('metadata/reduced_covariance_parallel_source_v2_fixture_restoration_review_20261003.json'),
        Path('metadata/full_uniform_covariance_qualification_plan_20261002.json'),
        Path('metadata/parallel_covariance_readback_plan_20261003_v1.json'),Path(dependency),Path('/usr/bin/prlimit'),Path(sys.executable)]
    pins={str(p):sha(p) for p in paths};pins.update(bindings)
    expected=dict(logical_cases=75188,model_setting_rows=622080,unique_cohorts=4340,unique_designs=130200,
        audit_rows=1302000,setting_audit_links=6220800)
    pp='metadata/full_reduced_covariance_qualification_plan_20261003_v2.json';output='results/phylogeny/full-exact-retained-covariance-qualification-20261003-v2'
    assert not Path(output).exists()
    scope=('Complete original75188cases/622080settings/4340cohorts/130200designs:1302000exact retained covariance audits/'
        '6220800setting links, both modes and allfive working trees. Explicit separately closed parallel readback of '
        'the same original full seven-kernel producer required; original serial reader path neither required,written '
        'nor replaced. Source plan/producer/alternate-reader/worker/resource identity, exactV2operator/cone proof and '
        'complete original journals checked. Frozen reduction/independentgesvd/SQL and unchanged principal Gram/error '
        'envelopes; all review/nonready/constant cases retained. V2 software byte-restoration finding preserved; stronger '
        'V3software admission rehashes the entire positive fixture graph. New immutable source adapter/output/controller '
        'namespace; original serial/qualification/timing queues untouched. Full source/artifact/two new original journals '
        'before numerical closure. No tolerance relaxation,removed exception kernel,optimizer,model acceptance,'
        'nonuniform qualification,accepted phylogeny/reconciliation/calibration or completed biological aims.')
    plan=dict(qualification_plan='metadata/full_uniform_covariance_qualification_plan_20261002.json',
        qualification_completion='metadata/full_uniform_covariance_parallel_completed_20261003_v1.json',
        parallel_readback_plan='metadata/parallel_covariance_readback_plan_20261003_v1.json',
        exact_plan=exact_plan,exact_completion=exact_completion,expected=expected,
        trees=['mafft_guide','pmsf_mafft_profile','pmsf_profile_mafft','pmsf_profile_profile','profile_guide'],
        output=output,pins=pins,resources=resources,dependencies=[dependency],
        completion='metadata/full_reduced_covariance_qualification_completed_20261003_v2.json',
        launch_inventory='metadata/full_reduced_covariance_qualification_launches_20261003_v2.json',scope=scope)
    create(pp,plan);prefix=['/usr/bin/prlimit','--as='+str(12*2**30),'--cpu=21600','--fsize='+str(2**30),'--']
    producer=launch('reduced-covariance-parallel-source-v2',prefix+[sys.executable,'scripts/run_full_reduced_covariance_qualification_v2.py','--plan',pp],
        [dependency],pp,cpus=2,memory=16)
    reader=launch('reduced-covariance-parallel-source-v2-readback',prefix+[sys.executable,'scripts/run_full_reduced_covariance_qualification_v2.py','--plan',pp,'--reader'],
        [producer],pp,cpus=2,memory=16)
    cp='metadata/full_reduced_covariance_qualification_completion_plan_20261003_v2.json'
    create(cp,dict(source_plan=pp,producer_receipt=output+'/receipt.json',independent_readback=output+'/readback.json',
        producer_status=PRODUCER_STATUS,reader_status=READER_STATUS,
        completed_status='complete_verified_full_exact_retained_uniform_covariance_qualification_v2',
        summary_fields=SUMMARY_FIELDS,launches=[producer,reader],pins={p:sha(p) for p in [pp,producer,reader,
            'scripts/close_full_triad_sequence_stage.py','scripts/record_completed_process_handoffs_v2.py']},
        output=plan['completion'],scope=scope))
    closer=launch('reduced-covariance-parallel-source-v2-closure',[sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',cp],
        [producer,reader],cp,cpus=2,memory=16)
    create(plan['launch_inventory'],dict(status='queued_complete_retained_covariance_after_alternate_source_closure',
        source_plan=pp,source_plan_sha256=sha(pp),launches=[producer,reader,closer],expected=expected,
        separately_replayed_real_certificates=8680,real_certificate_basis_counts=counts,
        original_serial_jobs_changed=False,original_reader_path_replaced=False,production_fitting_launched=False,
        all_eight_aims_incomplete=True,gpu=False,new_cost_usd=0))
    print(json.dumps(resources),flush=True)


if __name__=='__main__':main()
