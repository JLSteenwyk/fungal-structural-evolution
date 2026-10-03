#!/usr/bin/env python3
"""Checkpoint complete qualified-input runtime probes and retain every case."""
import argparse
from collections import Counter
import fcntl
import gzip
import json
import os
from pathlib import Path
import shutil
import time
import psutil
from full_retained_shared_entity_fit_sources import load,cohorts,operators_for
from full_covariance_qualification_sources import folded_operators
from full_retained_shared_entity_timing import groups,probe,estimate
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def sources(plan,path):
    fitpath=Path(plan['fit_plan']);assert sha(fitpath)==plan['fit_plan_sha256']
    fit=json.loads(fitpath.read_text());assert fit['independent_backend']=='component_spectral_v1'
    assert fit['launch_state']=='not_launched_or_queued' and not Path(fit['output']).exists()
    assert plan['scaled_variance_points']==[0.,1.]
    source,bindings=load(fit,fitpath)
    for p,d in plan['pins'].items():bind(bindings,p,d)
    bind(bindings,path);verify(bindings)
    return fit,source,bindings


def atomic(path,value):
    temporary=path.with_suffix(path.suffix+'.partial')
    with temporary.open('w') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(temporary,path)


def validate_probes(records,selected,counts,plan,fit):
    assert len(records)==len(selected) and {r['group_id'] for r in records}==set(selected)
    for record in records:
        expected=selected[record['group_id']]
        assert record['representative']==expected['identity'] and record['eligible_candidates']==counts[record['group_id']]
        assert record['selection_active_columns']==expected['matrix'].shape[1]
        assert record['selection_condition']==expected['rank'][1]
        assert record['variance_points']==plan['scaled_variance_points'] and record['scientific_eligibility'] is False
        assert record['status'] in ['timed_all_declared_points_agree_only','timing_group_requires_review','timing_group_construction_requires_review']
        if record['status']=='timing_group_construction_requires_review':
            assert record['error_type'] in ['ValueError','ArithmeticError','LinAlgError'];continue
        assert [p['scaled_variance'] for p in record['points']]==plan['scaled_variance_points']
        for name in ['source_validation_seconds','primary_constructor_seconds','production_qualification_guard_seconds','reader_qualification_seconds','independent_constructor_seconds']:
            assert isinstance(record[name],(int,float)) and 0<=record[name]<float('inf')
        assert record['parameter_names']==expected['source_audit']['retained_kernel_names'][1:]
        assert len(record['kernel_normalization'])==len(record['parameter_names'])
        assert all(isinstance(v,(int,float)) and 0<v<float('inf') for v in record['kernel_normalization'])
        for p in record['points']:
            assert p['status'] in ['timed_likelihood_agreement_only','timing_probe_numerical_agreement_requires_review','timing_probe_precision_requires_review']
            if p['status']=='timing_probe_precision_requires_review':
                assert p['error_type'] in ['ValueError','ArithmeticError','LinAlgError'];continue
            for name in ['primary_evaluation_seconds','independent_evaluation_seconds','objective_absolute_error','maximum_coordinate_gradient_error']:
                assert isinstance(p[name],(int,float)) and 0<=p[name]<float('inf')
            comparisons=['coefficient_and_conditional_covariance_comparison_passed','profiled_scale_comparison_passed','variance_components_comparison_passed']
            assert all(type(p[k]) is bool for k in comparisons)
            agreement=(p['objective_absolute_error']<=fit['independent_audit']['replay']['objective_atol'] and
                p['maximum_coordinate_gradient_error']<=fit['independent_audit']['replay']['gradient_atol'] and all(p[k] for k in comparisons))
            assert (p['status']=='timed_likelihood_agreement_only') is agreement
        assert (record['status']=='timed_all_declared_points_agree_only')==all(p['status']=='timed_likelihood_agreement_only' for p in record['points'])


def run(path,stop_after_cohorts=None):
    started=time.perf_counter();plan=json.loads(path.read_text());fit,source,bindings=sources(plan,path)
    root=Path(plan['output']);assert shutil.disk_usage(root.parent).free>=plan['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(exist_ok=True);lock=(root/'stage.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (root/'receipt.json').exists(),'Completed timing stage cannot restart'
    stage=dict(schema='full-retained-shared-entity-input-timing-v1',plan_sha256=sha(path),fit_contract=source['fit_contract'])
    marker=root/'stage_plan.json'
    if marker.exists():assert json.loads(marker.read_text())==stage
    else:atomic(marker,stage)
    folder=root/'cohorts';folder.mkdir(exist_ok=True);manifest=[];all_probes=[];statuses=Counter();total=eligible=0
    for number,(cohort,rows,entries) in enumerate(cohorts(source,fit),1):
        assert shutil.disk_usage(root).free>=plan['resources']['minimum_free_disk_gib']*2**30
        census,selected,counts=groups(source,fit,cohort,entries);cid=cohort['cohort_id']
        cp=folder/(cid+'.receipt.json');fp=folder/(cid+'.census.jsonl.gz');pp=folder/(cid+'.probes.json')
        if cp.exists():
            saved=json.loads(cp.read_text());assert saved['stage']==stage and saved['cohort_id']==cid
            assert saved['census_sha256']==sha(fp) and saved['probes_sha256']==sha(pp)
            with gzip.open(fp,'rt') as f:assert [json.loads(line) for line in f]==census
            probes=json.loads(pp.read_text());validate_probes(probes,selected,counts,plan,fit)
        else:
            operators={}
            for r in selected.values():
                mode=r['key'][0]
                if mode not in operators:operators[mode]=operators_for(source,rows,r['source_audit'])
            probes=[probe(source,fit,plan,rows,operators[r['key'][0]],r,counts[g],g) for g,r in sorted(selected.items())]
            validate_probes(probes,selected,counts,plan,fit)
            temp=fp.with_suffix('.gz.partial')
            with gzip.open(temp,'wt') as f:
                for r in census:f.write(json.dumps(r,sort_keys=True,allow_nan=False)+'\n')
            with temp.open('rb') as f:os.fsync(f.fileno())
            os.replace(temp,fp);atomic(pp,probes)
            saved=dict(stage=stage,cohort_id=cid,census_rows=len(census),eligible_candidates=sum(counts.values()),
                timing_groups=len(probes),census_sha256=sha(fp),probes_sha256=sha(pp))
            atomic(cp,saved)
        assert saved['census_rows']==len(census) and saved['eligible_candidates']==sum(counts.values()) and saved['timing_groups']==len(selected)
        manifest.append(dict(cohort_id=cid,census_path=str(fp.relative_to(root)),census_sha256=sha(fp),
            probes_path=str(pp.relative_to(root)),probes_sha256=sha(pp),receipt_path=str(cp.relative_to(root)),receipt_sha256=sha(cp)))
        total+=len(census);eligible+=sum(counts.values());statuses.update(r['source_disposition'] for r in census);all_probes.extend(probes)
        print('full_shared_entity_timing_cohort',number,'/',len(source['cohorts']),'cases',len(census),'timing_groups',len(probes),flush=True)
        if stop_after_cohorts==number:raise InterruptedError('Software timing checkpoint contract')
    expected=source['qualification_completion']['unique_designs']*2*len(fit['methods'])*len(fit['loading_modes'])*len(fit['trees'])
    assert total==expected
    planning=estimate(all_probes,fit);assert planning['measured_candidate_coverage']+planning['unmeasured_review_candidate_coverage']==eligible
    atomic(root/'cohort_manifest.json',manifest);atomic(root/'conditional_planning.json',planning)
    summary=dict(logical_cases=len(source['ids']),model_setting_rows=source['qualification_completion']['model_setting_rows'],
        unique_cohorts=len(manifest),candidate_rows=total,eligible_candidates=eligible,timing_groups=len(all_probes),
        source_status_counts=dict(statuses),timing_status_counts=dict(Counter(p['status'] for p in all_probes)),
        measured_candidate_coverage=planning['measured_candidate_coverage'],unmeasured_review_candidate_coverage=planning['unmeasured_review_candidate_coverage'],
        conditional_budget_weighted_seconds=planning['conditional_budget_weighted_seconds'])
    verify(bindings)
    artifacts={str(p.relative_to(root)):sha(p) for p in [marker,root/'cohort_manifest.json',root/'conditional_planning.json',
        *[root/m[k] for m in manifest for k in ['census_path','probes_path','receipt_path']]]}
    result=dict(status='complete_full_retained_shared_entity_timing_pending_accounting_readback_v1',plan_sha256=sha(path),
        fit_contract=source['fit_contract'],**summary,source_hashes=bindings,artifacts=artifacts,
        observed_wall_seconds=time.perf_counter()-started,observed_current_rss_bytes=psutil.Process().memory_info().rss,
        production_finish_eta=None,scientific_eligibility=False,scope=plan['scope'])
    atomic(root/'receipt.json',result);print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();run(a.plan)
