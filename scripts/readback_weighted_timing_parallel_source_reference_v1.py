#!/usr/bin/env python3
"""Read back every source case, deterministic probe choice and planning sum."""
import argparse
from collections import Counter
import fcntl
import gzip
import json
from pathlib import Path
from full_weighted_shared_entity_fit_sources_parallel_v1 import cohorts,operators_for
from full_weighted_shared_entity_timing import groups,estimate
from prepare_weighted_timing_parallel_source_reference_v1 import sources,SCHEMA,PRODUCER,READER
from full_weighted_timing_contracts_v2 import validate_probes,replay,review_inputs
from full_weighted_fit_exports import atomic
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def run(path,output):
    assert not Path(output).exists(), 'Completed readback cannot restart'
    root=Path(json.loads(path.read_text())['output'])
    lock=(root/'reader.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (root/'reader_completed.json').exists(), 'Completed reader cannot restart with another output'
    plan=json.loads(path.read_text());fit,source,bindings=sources(plan,path);original=dict(bindings)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text());bind(bindings,rp)
    assert receipt['status']==PRODUCER and receipt['plan_sha256']==sha(path)
    assert receipt['source_hashes']==original and receipt['fit_contract']==source['fit_contract'] and receipt['scientific_eligibility'] is False
    for name,d in receipt['artifacts'].items():bind(bindings,root/name,d)
    verify(bindings)
    stage=dict(schema=SCHEMA,plan_sha256=sha(path),fit_contract=source['fit_contract'])
    assert json.loads((root/'stage_plan.json').read_text())==stage
    manifest=json.loads((root/'cohort_manifest.json').read_text());assert [m['cohort_id'] for m in manifest]==[c['cohort_id'] for c in source['cohorts']]
    all_probes=[];counts=Counter();total=eligible=0;reviews=[]
    for (cohort,rows,entries),part in zip(cohorts(source,fit),manifest):
        census,selected,multiplicity=groups(source,fit,cohort,entries)
        for key in ['census','probes','receipt']:assert sha(root/part[key+'_path'])==part[key+'_sha256']
        with gzip.open(root/part['census_path'],'rt') as f:assert [json.loads(line) for line in f]==census
        probes=json.loads((root/part['probes_path']).read_text());validate_probes(probes,selected,multiplicity,plan,fit)
        replay(source,fit,plan,rows,probes,selected,multiplicity)
        for p in probes:reviews.extend(review_inputs(root,p,selected[p['group_id']],source,rows,read=True))
        checkpoint=json.loads((root/part['receipt_path']).read_text())
        assert checkpoint==dict(stage=stage,cohort_id=cohort['cohort_id'],census_rows=len(census),eligible_candidates=sum(multiplicity.values()),
            timing_groups=len(probes),census_sha256=part['census_sha256'],probes_sha256=part['probes_sha256'])
        total+=len(census);eligible+=sum(multiplicity.values());counts.update(c['source_disposition'] for c in census);all_probes.extend(probes)
    assert total==source['numerical_completion']['designs']*2*len(fit['methods'])*len(fit['loading_modes'])*len(fit['trees'])*len(fit['policies'])==fit['expected']['candidate_rows']
    planning=estimate(all_probes,fit);assert planning==json.loads((root/'conditional_planning.json').read_text())
    summary=dict(logical_cases=len(source['ids']),model_setting_rows=source['numerical_completion']['settings'],
        unique_cohorts=len(manifest),candidate_rows=total,eligible_candidates=eligible,timing_groups=len(all_probes),
        source_status_counts=dict(counts),timing_status_counts=dict(Counter(p['status'] for p in all_probes)),
        measured_candidate_coverage=planning['measured_candidate_coverage'],unmeasured_review_candidate_coverage=planning['unmeasured_review_candidate_coverage'],
        conditional_budget_weighted_seconds=planning['conditional_budget_weighted_seconds'])
    assert all(receipt[k]==v for k,v in summary.items()) and receipt['production_finish_eta'] is None
    paths=[root/'stage_plan.json',root/'cohort_manifest.json',root/'conditional_planning.json',
        *[root/m[k] for m in manifest for k in ['census_path','probes_path','receipt_path']],*reviews]
    assert receipt['artifacts']=={str(p.relative_to(root)):sha(p) for p in paths}
    assert receipt['fits_computed']==0 and receipt['nonuniform_weighting_accepted'] is False and receipt['component_variance_attribution_accepted'] is False
    verify(bindings)
    result=dict(status=READER,plan_sha256=sha(path),
        producer_receipt_sha256=sha(rp),fit_contract=source['fit_contract'],**summary,source_hashes=bindings,
        fresh_numeric_probes_replayed=len(all_probes),independent_hardware_timing_reimplementation=False,fits_computed=0,
        nonuniform_weighting_accepted=False,component_variance_attribution_accepted=False,production_finish_eta=None,scientific_eligibility=False,scope=plan['scope'])
    atomic(Path(output),result)
    atomic(root/'reader_completed.json',dict(output=str(output),sha256=sha(output),status=READER,plan_sha256=sha(path)))
    print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.plan,a.output)
