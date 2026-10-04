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
from full_weighted_shared_entity_fit_sources import load,cohorts,operators_for
from full_weighted_shared_entity_timing import groups,probe,estimate,SCHEMA,PRODUCER
from full_weighted_timing_contracts import validate_probes,replay,review_inputs
from full_weighted_fit_exports import atomic,finish_gzip,failure_capture
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


def run(path,stop_after_cohorts=None):
    started=time.perf_counter();plan=json.loads(path.read_text());fit,source,bindings=sources(plan,path)
    root=Path(plan['output']);assert shutil.disk_usage(root.parent).free>=plan['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(exist_ok=True);lock=(root/'stage.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (root/'receipt.json').exists(),'Completed timing stage cannot restart'
    assert not (root/'failures').exists(), 'Failed timing stage preserved; use a fresh version'
    stage=dict(schema=SCHEMA,plan_sha256=sha(path),fit_contract=source['fit_contract'])
    marker=root/'stage_plan.json'
    if marker.exists():assert json.loads(marker.read_text())==stage
    else:atomic(marker,stage)
    folder=root/'cohorts';folder.mkdir(exist_ok=True);manifest=[];all_probes=[];statuses=Counter();total=eligible=0;reviews=[]
    for number,(cohort,rows,entries) in enumerate(cohorts(source,fit),1):
        assert shutil.disk_usage(root).free>=plan['resources']['minimum_free_disk_gib']*2**30
        census,selected,counts=groups(source,fit,cohort,entries);cid=cohort['cohort_id']
        cp=folder/(cid+'.receipt.json');fp=folder/(cid+'.census.jsonl.gz');pp=folder/(cid+'.probes.json')
        if cp.exists():
            saved=json.loads(cp.read_text());assert saved['stage']==stage and saved['cohort_id']==cid
            assert saved['census_sha256']==sha(fp) and saved['probes_sha256']==sha(pp)
            with gzip.open(fp,'rt') as f:assert [json.loads(line) for line in f]==census
            probes=json.loads(pp.read_text());validate_probes(probes,selected,counts,plan,fit)
            replay(source,fit,plan,rows,probes,selected,counts)
        else:
            assert not fp.exists() and not pp.exists(), 'Unclosed artifacts preserved; use a fresh stage version'
            operators={};probes=[]
            for g,r in sorted(selected.items()):
                key=(r['identity']['loading_mode'],tuple(r['source_audit']['retained_kernel_names']))
                if key not in operators:operators[key]=operators_for(source,rows,r['source_audit'])
                try:probes.append(probe(source,fit,plan,rows,operators[key],r,counts[g],g))
                except Exception as error:
                    failure_capture(root/'failures',r['identity'],probes,rows,r['diagonal'],r['matrix'],r['response'],
                        source['factors'][r['identity']['tree']][rows],error)
                    raise
            validate_probes(probes,selected,counts,plan,fit)
            temp=fp.with_suffix('.gz.partial')
            with gzip.open(temp,'xt') as f:
                for r in census:f.write(json.dumps(r,sort_keys=True,allow_nan=False)+'\n')
            with temp.open('rb') as f:os.fsync(f.fileno())
            finish_gzip(temp,fp);atomic(pp,probes)
            for record in probes:review_inputs(root,record,selected[record['group_id']],source,rows)
            saved=dict(stage=stage,cohort_id=cid,census_rows=len(census),eligible_candidates=sum(counts.values()),
                timing_groups=len(probes),census_sha256=sha(fp),probes_sha256=sha(pp))
            atomic(cp,saved)
        for record in probes:reviews.extend(review_inputs(root,record,selected[record['group_id']],source,rows,read=True))
        assert saved['census_rows']==len(census) and saved['eligible_candidates']==sum(counts.values()) and saved['timing_groups']==len(selected)
        manifest.append(dict(cohort_id=cid,census_path=str(fp.relative_to(root)),census_sha256=sha(fp),
            probes_path=str(pp.relative_to(root)),probes_sha256=sha(pp),receipt_path=str(cp.relative_to(root)),receipt_sha256=sha(cp)))
        total+=len(census);eligible+=sum(counts.values());statuses.update(r['source_disposition'] for r in census);all_probes.extend(probes)
        print('full_shared_entity_timing_cohort',number,'/',len(source['cohorts']),'cases',len(census),'timing_groups',len(probes),flush=True)
        if stop_after_cohorts==number:raise InterruptedError('Software timing checkpoint contract')
    expected=source['numerical_completion']['designs']*2*len(fit['methods'])*len(fit['loading_modes'])*len(fit['trees'])*len(fit['policies'])
    assert total==expected==fit['expected']['candidate_rows']
    planning=estimate(all_probes,fit);assert planning['measured_candidate_coverage']+planning['unmeasured_review_candidate_coverage']==eligible
    atomic(root/'cohort_manifest.json',manifest);atomic(root/'conditional_planning.json',planning)
    summary=dict(logical_cases=len(source['ids']),model_setting_rows=source['numerical_completion']['settings'],
        unique_cohorts=len(manifest),candidate_rows=total,eligible_candidates=eligible,timing_groups=len(all_probes),
        source_status_counts=dict(statuses),timing_status_counts=dict(Counter(p['status'] for p in all_probes)),
        measured_candidate_coverage=planning['measured_candidate_coverage'],unmeasured_review_candidate_coverage=planning['unmeasured_review_candidate_coverage'],
        conditional_budget_weighted_seconds=planning['conditional_budget_weighted_seconds'])
    verify(bindings)
    artifacts={str(p.relative_to(root)):sha(p) for p in [marker,root/'cohort_manifest.json',root/'conditional_planning.json',
        *[root/m[k] for m in manifest for k in ['census_path','probes_path','receipt_path']],*reviews]}
    result=dict(status=PRODUCER,plan_sha256=sha(path),
        fit_contract=source['fit_contract'],**summary,source_hashes=bindings,artifacts=artifacts,
        observed_wall_seconds=time.perf_counter()-started,observed_current_rss_bytes=psutil.Process().memory_info().rss,
        production_finish_eta=None,fits_computed=0,nonuniform_weighting_accepted=False,component_variance_attribution_accepted=False,scientific_eligibility=False,scope=plan['scope'])
    atomic(root/'receipt.json',result);print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();run(a.plan)
