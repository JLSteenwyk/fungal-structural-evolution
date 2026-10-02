#!/usr/bin/env python3
"""Integrate every closed optimization follow-up and export the complete matched grid."""
import argparse
from collections import Counter
from contextlib import contextmanager
import fcntl
import gzip
import io
import json
import math
from pathlib import Path
import shutil
from full_whole_protein_comparison_sources import load_sources
from export_whole_protein_ml_comparisons import comparison,PASS,REFINED
from run_whole_protein_ml import digest
from reference_measurement_union_sources import bind,verify
from screen_duplication_alignment_reuse import sha


@contextmanager
def deterministic_gzip(path):
    # Parameter hashes participate in restart checkpoints: omit time and filename.
    with Path(path).open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',compresslevel=1,fileobj=raw,mtime=0) as compressed:
            with io.TextIOWrapper(compressed,encoding='utf-8') as output:
                yield output


def selected(entry,saved,follow,follow_root,support):
    key=entry['fit_input_id'],entry['tree'];spec=saved['specification'];original=saved['payload'];status=original['status']
    assert (saved['fit_input_id'],saved['tree'])==key and status==entry['status'] and saved['payload_sha256']==digest(original)
    chosen=original;kind='original_production';source=entry['path'];source_hash=entry['sha256'];effective=status;case_status=None;case_path=None;case_hash=None
    theta=original.get('log1p_ratios');raw=original
    if status=='ml_candidate_requires_optimization_review':
        row=follow[key];case_path=str(Path(follow_root)/row['path']);case_hash=row['sha256'];assert sha(case_path)==case_hash
        case=json.loads(Path(case_path).read_text());identity=case['identity'];result=case['result'];case_status=result['status']
        assert identity['fit_input_id']==key[0] and identity['tree']==key[1] and identity['original_fit']==entry['path'] and identity['original_fit_sha256']==entry['sha256'] and identity['input_sha256']==saved['input_sha256']
        assert case['result_sha256']==digest(result) and case['original_status']==status and case['scientific_eligibility'] is False and case_status==row['status']
        if case_status=='numerical_followup_error_requires_review':kind='original_retained_after_followup_error'
        else:
            assert case_status in ['numerical_followup_passed_pending_independent_readback','numerical_followup_requires_review']
            chosen=result['fitted'];raw=result;theta=result['selected_theta'];source=case_path;source_hash=case_hash
            kind='full_grid_refinement' if result['recovery'] is None else 'full_grid_recovery'
            effective=REFINED if case_status=='numerical_followup_passed_pending_independent_readback' else 'ml_full_grid_followup_requires_review'
    else:assert key not in follow and status in [PASS,'fit_error_requires_review']
    error=status=='fit_error_requires_review';p=len(spec['columns'])-1
    if not error:
        assert math.isfinite(chosen['negative_profiled_ml']) and len(chosen['beta'])==p and chosen['variance_profile_denominator']==spec['records']
    supported=support.resolve(key[0],saved['input_sha256'])
    summary=dict(fit_input_id=key[0],tree=key[1],input_sha256=saved['input_sha256'],original_status=status,effective_status=effective,
        original_negative_profiled_ml=original.get('negative_profiled_ml'),negative_profiled_ml=None if error else chosen['negative_profiled_ml'],records=spec['records'],fixed_coefficients=p,
        original_production_fit=entry['path'],original_production_sha256=entry['sha256'],selected_source_kind=kind,selected_source_path=source,selected_source_sha256=source_hash,
        followup_status=case_status,followup_source_path=case_path,followup_source_sha256=case_hash,numerical_fit_verified=not error,
        support_geometry_id=supported['geometry_id'],support_classification=supported['classification'],original_support_classification=supported['original_classification'],support_source_kind=supported['source_kind'])
    params=dict(fit_input_id=key[0],tree=key[1],columns=spec['columns'],specification=spec,selected_source_path=source,selected_source_sha256=source_hash,effective_status=effective,numerical_fit_verified=not error,
        covariate_scales=None if error else original['covariate_scales'],scaled_beta=None if error else chosen['beta'],
        scaled_conditional_beta_covariance=None if error else chosen['conditional_beta_covariance'],raw_unit_beta=None if error else raw['raw_unit_beta'],
        raw_unit_conditional_beta_covariance=None if error else raw['raw_unit_conditional_beta_covariance'],selected_log1p_ratios=None if error else theta,
        profiled_scale=None if error else chosen['profiled_scale'],variance_profile_denominator=None if error else chosen['variance_profile_denominator'],scientific_eligibility=False)
    return summary,params


def run(path):
    plan=json.loads(Path(path).read_text())
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    original,follow,follow_root,lr,support,bindings=load_sources(plan,path)
    out=Path(plan['output']);out.mkdir(exist_ok=True,parents=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (out/'receipt.json').exists()
        config=dict(plan_sha256=sha(path),source_hashes=bindings);cp=out/'configuration.json'
        if cp.exists():assert json.loads(cp.read_text())==config
        else:
            temp=out/'configuration.tmp';temp.write_text(json.dumps(config,indent=2)+'\n');temp.replace(cp)
        summaries={};old_counts=Counter();counts=Counter();cases=Counter();kinds=Counter()
        sp=out/'fit_summaries.jsonl';param=out/'selected_fit_parameters.jsonl.gz'
        with sp.with_suffix('.tmp').open('w') as sf,deterministic_gzip(param.with_suffix('.tmp')) as pf:
            for key,entry in original.items():
                saved=json.loads(Path(entry['path']).read_text());s,p=selected(entry,saved,follow,follow_root,support);summaries[key]=s
                old_counts[s['original_status']]+=1;counts[s['effective_status']]+=1;kinds[s['selected_source_kind']]+=1
                if s['followup_status'] is not None:cases[s['followup_status']]+=1
                sf.write(json.dumps(s,separators=(',',':'),allow_nan=False)+'\n');pf.write(json.dumps(p,separators=(',',':'),allow_nan=False)+'\n')
                if len(summaries)%25000==0:print('Full optimized model summaries',len(summaries),flush=True)
        sp.with_suffix('.tmp').replace(sp);param.with_suffix('.tmp').replace(param)
        trees=sorted({k[1] for k in summaries});assert len(trees)==plan['expected']['trees']
        ids={k[0] for k in summaries};assert len(ids)==plan['expected']['unique_inputs'] and len(summaries)==plan['expected']['fit_summaries']
        assert all((i,t) in summaries for i in ids for t in trees) and sum(cases.values())==len(follow)
        artifacts={p.name:sha(p) for p in [sp,param]};counts_comparison=Counter();relations=Counter();total=0;reused=0
        for tree in trees:
            target=out/(tree+'.jsonl.gz');checkpoint=out/(tree+'.checkpoint.json');configuration=dict(plan_sha256=sha(path),summary_sha256=sha(sp),parameters_sha256=sha(param))
            if checkpoint.exists():
                c=json.loads(checkpoint.read_text());assert c['configuration']==configuration and sha(target)==c['sha256'];reused+=1
            else:
                n=0;local=Counter();rel=Counter();tmp=target.with_suffix('.tmp')
                with gzip.open(tmp,'wt',compresslevel=1) as f,(Path(plan['links'])/'comparison_input_map.jsonl').open() as source:
                    for line in source:
                        link=json.loads(line);left=summaries[link['left_input'],tree];right=summaries[link['right_input'],tree]
                        row=dict(**link,tree=tree,left=left,right=right,**comparison(left,right,link['relation']))
                        f.write(json.dumps(row,separators=(',',':'),allow_nan=False)+'\n');n+=1;local[row['comparison_status']]+=1;rel[link['relation']]+=1
                assert n==plan['expected']['comparisons_per_tree'] and dict(rel)==lr['relation_counts'];tmp.replace(target)
                c=dict(configuration=configuration,sha256=sha(target),comparisons=n,status_counts=dict(local),relation_counts=dict(rel));tmp=checkpoint.with_suffix('.tmp');tmp.write_text(json.dumps(c)+'\n');tmp.replace(checkpoint)
            assert c['comparisons']==plan['expected']['comparisons_per_tree'];total+=c['comparisons'];counts_comparison.update(c['status_counts']);relations.update(c['relation_counts'])
            artifacts[target.name]=sha(target);artifacts[checkpoint.name]=sha(checkpoint)
            assert sum(p.stat().st_size for p in out.iterdir() if p.is_file())<=plan['resources']['maximum_output_gib']*2**30 and shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
            print('Full optimized comparison tree',tree,total,flush=True)
        assert total==plan['expected']['comparisons'];bind(bindings,cp);verify(bindings)
        result=dict(status='complete_full_whole_protein_comparisons_pending_independent_readback',fit_summaries=len(summaries),parameter_records=len(summaries),comparisons=total,trees=trees,
            original_fit_status_counts=dict(old_counts),effective_fit_status_counts=dict(counts),followup_status_counts=dict(cases),selected_source_kind_counts=dict(kinds),
            comparison_status_counts=dict(counts_comparison),relation_counts=dict(relations),followup_cases_integrated=len(follow),checked_reused_tree_checkpoints=reused,
            plan_sha256=sha(path),source_hashes=bindings,artifacts=artifacts,scientific_eligibility=False,scope=plan['scope'])
        temp=out/'receipt.tmp';temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(out/'receipt.json')
        print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
