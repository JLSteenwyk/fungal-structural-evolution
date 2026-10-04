"""Admit unchanged fit inputs after full parallel timing/checkpoint closure.

This is a complete source/census/accounting gate. It does not repeat the closed
reader's numeric probes or interpret timings as inferential acceptance.
Operational configuration stays separate from the immutable measured fit plan.
"""
from collections import Counter
import gzip
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from full_exact_covariance_sources import closed_subset
from full_weighted_shared_entity_fit_sources_parallel_v1 import load, cohorts
from full_weighted_shared_entity_timing import groups, estimate
from full_weighted_timing_contracts_v2 import validate_probes, review_inputs
from prepare_weighted_timing_parallel_source_reference_v1 import SCHEMA, PRODUCER, READER
from reference_measurement_union_sources import bind, verify

COMPLETED='complete_verified_full_four_control_input_timing_accounting_v2'
FIELDS=['logical_cases','model_setting_rows','unique_cohorts','candidate_rows','eligible_candidates',
    'timing_groups','source_status_counts','timing_status_counts','measured_candidate_coverage',
    'unmeasured_review_candidate_coverage','conditional_budget_weighted_seconds']
FULL=dict(logical_cases=75188,cohorts=4340,designs=130200,fit_inputs=260400,settings=622080,
    candidate_rows=20832000,setting_fit_links=49766400)


def checkpoint_pair(producer,reader,part,index,census,probes,multiplicity,planning,reviews):
    """Verify both complete guards/counts/identities; this is not native replay."""
    expected=dict(cohort_index=index,cohort_id=part['cohort_id'],manifest=part,
        candidate_rows=len(census),eligible_candidates=sum(multiplicity.values()),
        timing_groups=len(probes),source_status_counts=dict(Counter(r['source_disposition'] for r in census)),
        timing_status_counts=dict(Counter(r['status'] for r in probes)),planning=planning,review_artifacts=reviews)
    for record in [producer,reader]:
        assert all(record[k]==v for k,v in expected.items())
        assert record['cached_numeric_inputs_preserved'] is True
        assert record['submitted_numeric_inputs_preserved'] is True
        assert record['scientific_eligibility'] is False
        assert record['worker_address_space_limit_bytes']==12*2**30
        for key in ['cached_input_array_bindings_checked','submitted_task_array_bindings_checked']:
            assert type(record[key]) is int and record[key]>0
        worker=record['worker']
        assert type(worker['pid']) is int and worker['pid']>0
        assert type(worker['created']) in [int,float] and worker['created']>0
        assert isinstance(worker['cmdline'],list) and worker['cmdline'] and all(type(v) is str for v in worker['cmdline'])
    assert producer['review_artifacts']==reader['review_artifacts']


def closed_timing(fit_path, timing_path, completion_path):
    fit_path,timing_path,completion_path=map(Path,[fit_path,timing_path,completion_path])
    fit=json.loads(fit_path.read_text());timing=json.loads(timing_path.read_text())
    assert timing['fit_plan']==str(fit_path) and timing['fit_plan_sha256']==sha(fit_path)
    if 'completion' in timing:assert timing['completion']==str(completion_path)
    assert timing['scaled_variance_points']==[0.,1.] and fit['independent_backend']=='component_spectral_v1'
    assert fit['launch_state']=='not_launched_or_queued', 'Measured historical configuration must remain unchanged'
    assert completion_path.is_file(), 'Complete original timing closure is required before fitting admission'
    assert Path(fit['qualification_completion']).is_file(), 'Complete original weighted numerical closure is required'
    source,original=load(fit,fit_path)
    expected_sources=dict(original)
    for p,d in timing['pins'].items():bind(expected_sources,p,d)
    bind(expected_sources,timing_path);verify(expected_sources)
    root=Path(timing['output']);rp=root/'receipt.json';rb=root/'readback.json'
    producer=json.loads(rp.read_text());reader=json.loads(rb.read_text())
    assert producer['status']==PRODUCER and reader['status']==READER
    assert producer['plan_sha256']==reader['plan_sha256']==sha(timing_path)
    assert producer['fit_contract']==reader['fit_contract']==source['fit_contract']
    assert producer['source_hashes']==expected_sources
    assert reader['producer_receipt_sha256']==sha(rp)
    for r in [producer,reader]:
        assert r['scientific_eligibility'] is r['nonuniform_weighting_accepted'] is r['component_variance_attribution_accepted'] is False
        assert r['fits_computed']==0 and r['production_finish_eta'] is None
    assert reader['independent_hardware_timing_reimplementation'] is False
    assert reader['fresh_numeric_probes_replayed']==reader['timing_groups']
    bound_reader=dict(expected_sources);bind(bound_reader,rp)
    for name,d in producer['artifacts'].items():
        assert Path(name).is_relative_to('.') and not Path(name).is_absolute() and '..' not in Path(name).parts
        bind(bound_reader,root/name,d)
    assert reader['source_hashes']==bound_reader;verify(bound_reader)
    reader_checkpoints=reader['reader_checkpoint_artifacts']
    assert set(reader_checkpoints)=={'checkpoints/'+str(i).zfill(5)+'.reader.json' for i in range(producer['unique_cohorts'])}
    bindings=dict(bound_reader);bind(bindings,rb)
    for name,digest in reader_checkpoints.items():bind(bindings,root/name,digest)
    completion=closed_subset(completion_path,COMPLETED,
        [timing_path,rp,rb,*[root/name for name in producer['artifacts']],*[root/name for name in reader_checkpoints]],bindings)
    assert completion['producer_receipt']==str(rp) and completion['producer_receipt_sha256']==sha(rp)
    assert completion['independent_readback']==str(rb) and completion['independent_readback_sha256']==sha(rb)
    assert all(producer[k]==reader[k]==completion[k] for k in FIELDS)
    stage=dict(schema=SCHEMA,plan_sha256=sha(timing_path),fit_contract=source['fit_contract'])
    assert json.loads((root/'stage_plan.json').read_text())==stage
    manifest=json.loads((root/'cohort_manifest.json').read_text())
    assert [r['cohort_id'] for r in manifest]==[c['cohort_id'] for c in source['cohorts']]
    planning=estimate([],fit);statuses=Counter();timing_statuses=Counter();count=eligible=group_count=0
    paths=[root/'stage_plan.json',root/'cohort_manifest.json',root/'conditional_planning.json']
    for index,((cohort,rows,entries),part) in enumerate(zip(cohorts(source,fit),manifest)):
        census,selected,multiplicity=groups(source,fit,cohort,entries)
        for key in ['census','probes','receipt']:
            assert part[key+'_path']=='cohorts/'+cohort['cohort_id']+'.'+{'census':'census.jsonl.gz','probes':'probes.json','receipt':'receipt.json'}[key]
            assert sha(root/part[key+'_path'])==part[key+'_sha256']
            paths.append(root/part[key+'_path'])
        with gzip.open(root/part['census_path'],'rt') as f:assert [json.loads(l) for l in f]==census
        probes=json.loads((root/part['probes_path']).read_text());validate_probes(probes,selected,multiplicity,timing,fit)
        reviews=[]
        for record in probes:reviews.extend(review_inputs(root,record,selected[record['group_id']],source,rows,read=True))
        paths.extend(reviews)
        assert json.loads((root/part['receipt_path']).read_text())==dict(stage=stage,cohort_id=cohort['cohort_id'],
            census_rows=len(census),eligible_candidates=sum(multiplicity.values()),timing_groups=len(probes),
            census_sha256=part['census_sha256'],probes_sha256=part['probes_sha256'])
        local=estimate(probes,fit)
        wp=root/'checkpoints'/(str(index).zfill(5)+'.producer.json')
        wr=root/'checkpoints'/(str(index).zfill(5)+'.reader.json')
        checkpoint_pair(json.loads(wp.read_text()),json.loads(wr.read_text()),part,index,census,probes,multiplicity,local,
            {str(path.relative_to(root)):sha(path) for path in reviews})
        assert producer['artifacts'][str(wp.relative_to(root))]==sha(wp)
        assert reader_checkpoints[str(wr.relative_to(root))]==sha(wr)
        paths.append(wp)
        for field in ['measured_candidate_coverage','unmeasured_review_candidate_coverage']:planning[field]+=local[field]
        planning['group_planning_costs'].extend(local['group_planning_costs'])
        count+=len(census);eligible+=sum(multiplicity.values());group_count+=len(probes)
        statuses.update(r['source_disposition'] for r in census);timing_statuses.update(r['status'] for r in probes)
    assert count==fit['expected']['candidate_rows']==source['numerical_completion']['designs']*2*4*2*5*2
    assert producer['artifacts']=={str(p.relative_to(root)):sha(p) for p in paths}
    planning['conditional_budget_weighted_seconds']=sum(p['conditional_budget_weighted_producer_seconds']+
        p['conditional_budget_weighted_reader_seconds'] for p in planning['group_planning_costs'])
    assert planning==json.loads((root/'conditional_planning.json').read_text())
    assert planning['measured_candidate_coverage']+planning['unmeasured_review_candidate_coverage']==eligible
    summary=dict(logical_cases=len(source['ids']),model_setting_rows=source['numerical_completion']['settings'],
        unique_cohorts=len(manifest),candidate_rows=count,eligible_candidates=eligible,timing_groups=group_count,
        source_status_counts=dict(statuses),timing_status_counts=dict(timing_statuses),
        measured_candidate_coverage=planning['measured_candidate_coverage'],
        unmeasured_review_candidate_coverage=planning['unmeasured_review_candidate_coverage'],
        conditional_budget_weighted_seconds=planning['conditional_budget_weighted_seconds'])
    assert all(reader[k]==v for k,v in summary.items())
    marker=root/'reader_completed.json'
    assert json.loads(marker.read_text())==dict(output=str(rb),sha256=sha(rb),status=READER,plan_sha256=sha(timing_path))
    bind(bindings,marker);verify(bindings)
    return source,bindings,dict(**summary,fit_contract=source['fit_contract'],planning=planning,
        original_fit_plan=str(fit_path),original_fit_plan_sha256=sha(fit_path),timing_plan=str(timing_path),
        timing_completion=str(completion_path),timing_completion_sha256=sha(completion_path),
        fresh_admission_candidate_rows=count,closed_numeric_probe_groups=reader['fresh_numeric_probes_replayed'],
        admission_repeats_numeric_probes=False,full_parallel_checkpoint_pairs_checked=len(manifest),
        admission_retains_all_probe_arrays_in_memory=False,production_finish_eta=None,scientific_eligibility=False)


def production_contract(fit,admission):
    """Synthetic gate tests cannot become a reduced production launch scope."""
    assert fit['expected']==FULL
    assert (admission['logical_cases'],admission['unique_cohorts'],admission['candidate_rows'],admission['model_setting_rows'])==(75188,4340,20832000,622080)
    assert fit['methods']==['ml','reml'] and fit['loading_modes']==['signed','unsigned']
    assert fit['policies']==['uniform','background_node','background_pair','family_component'] and len(fit['trees'])==5
    assert fit['optimizer']==dict(gradient_tolerance=1e-6,max_evaluations=500,max_iterations=200,maximum_scaled_variance=1e6)
    assert fit['independent_audit']['column_batch']==32
    assert fit['independent_audit']['replay']==dict(gradient_atol=1e-6,objective_atol=1e-7)
    assert fit['independent_audit']['curvature']==dict(coordinate_step=.0001,gradient_atol=1e-6)
    assert fit['independent_audit']['optimizer']==dict(gradient_tolerance=1e-6,max_evaluations=500,max_iterations=200,objective_tolerance=1e-7)


def capacity(admission,fit):
    producer=sum(p['conditional_budget_weighted_producer_seconds'] for p in admission['planning']['group_planning_costs'])
    reader=sum(p['conditional_budget_weighted_reader_seconds'] for p in admission['planning']['group_planning_costs'])
    import math
    import numpy as np
    assert producer>=0 and reader>=0 and math.isclose(producer+reader,admission['conditional_budget_weighted_seconds'],
        rel_tol=8*np.finfo(float).eps*max(1,len(admission['planning']['group_planning_costs'])),abs_tol=1e-12)
    # Capacity only. Include an explicit uncalibrated allowance for source,
    # export and review paths, without claiming the point observations bound
    # any optimizer trajectory. Every source candidate still reaches exports.
    resources=dict(cpus=2,memory_gib=32,swap_gib=0,blas_threads=1,address_space_gib=24,
        per_file_limit_gib=2,minimum_available_ram_gib=64,minimum_free_disk_gib=fit['resources']['minimum_free_disk_gib'],
        output_scratch_reserve_gib=fit['resources']['output_scratch_reserve_gib'],
        conditional_budget_weighted_producer_seconds=producer,conditional_budget_weighted_reader_seconds=reader,
        producer_native_cpu_allocation_seconds=math.ceil(4*producer+604800),reader_native_cpu_allocation_seconds=math.ceil(4*reader+604800),
        source_export_review_capacity_margin_seconds=604800,observed_cost_capacity_multiplier=4,
        unmeasured_review_candidate_coverage=admission['unmeasured_review_candidate_coverage'],
        mathematical_runtime_bound=False,complete_runtime_calibration=False,production_finish_eta=None,
        wall_time_limit_enforced=False,new_cost_usd=0,gpu=False,resources_installed=False,
        scope='Conditional actual point costs receive fourfold capacity allowance plus168CPUhours per-role for uncalibrated source/export/review work. This is a resource allocation, not expected runtime or a convergence/finish guarantee. Unmeasured precision reviews are retained; all original candidates still require full exports and independent numeric/link readback. Local two-CPU/32GiB/no-swap/one-BLAS resources with separate producer/reader CPU allocations; no GPU or new charge.')
    assert resources['producer_native_cpu_allocation_seconds']<2**63 and resources['reader_native_cpu_allocation_seconds']<2**63
    return resources
