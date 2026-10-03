#!/usr/bin/env python3
"""Read back every source case, deterministic probe choice and planning sum."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
from full_shared_entity_fit_sources import cohorts
from full_shared_entity_timing import groups,estimate,probe
from full_covariance_qualification_sources import folded_operators
from prepare_full_shared_entity_timing import sources,validate_probes
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def run(path,output):
    plan=json.loads(path.read_text());fit,source,bindings=sources(plan,path);original=dict(bindings)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text());bind(bindings,rp)
    assert receipt['status']=='complete_full_shared_entity_timing_pending_accounting_readback' and receipt['plan_sha256']==sha(path)
    assert receipt['source_hashes']==original and receipt['fit_contract']==source['fit_contract'] and receipt['scientific_eligibility'] is False
    for name,d in receipt['artifacts'].items():bind(bindings,root/name,d)
    verify(bindings)
    stage=dict(schema='full-shared-entity-input-timing-v1',plan_sha256=sha(path),fit_contract=source['fit_contract'])
    assert json.loads((root/'stage_plan.json').read_text())==stage
    manifest=json.loads((root/'cohort_manifest.json').read_text());assert [m['cohort_id'] for m in manifest]==[c['cohort_id'] for c in source['cohorts']]
    all_probes=[];counts=Counter();total=eligible=0
    for (cohort,rows,entries),part in zip(cohorts(source,fit),manifest):
        census,selected,multiplicity=groups(source,fit,cohort,entries)
        for key in ['census','probes','receipt']:assert sha(root/part[key+'_path'])==part[key+'_sha256']
        with gzip.open(root/part['census_path'],'rt') as f:assert [json.loads(line) for line in f]==census
        probes=json.loads((root/part['probes_path']).read_text());validate_probes(probes,selected,multiplicity,plan,fit)
        # Reproduce every claimed precision/construction review instead of
        # permitting an invented failure to conceal an eligible timing group.
        operators={}
        for p in probes:
            if p['status']=='timed_all_declared_points_agree_only':continue
            representative=selected[p['group_id']];mode=representative['key'][0]
            if mode not in operators:operators[mode]=folded_operators(source,rows,mode)
            repeated=probe(source,fit,plan,rows,operators[mode],representative,multiplicity[p['group_id']],p['group_id'])
            assert repeated['status']==p['status']
            if p['status']=='timing_group_construction_requires_review':
                assert repeated['error_type']==p['error_type'] and repeated['error_message']==p['error_message']
            else:
                assert [v['status'] for v in repeated['points']]==[v['status'] for v in p['points']]
                for a,b in zip(repeated['points'],p['points']):
                    if a['status']=='timing_probe_precision_requires_review':assert a['error_type']==b['error_type'] and a['error_message']==b['error_message']
        checkpoint=json.loads((root/part['receipt_path']).read_text())
        assert checkpoint==dict(stage=stage,cohort_id=cohort['cohort_id'],census_rows=len(census),eligible_candidates=sum(multiplicity.values()),
            timing_groups=len(probes),census_sha256=part['census_sha256'],probes_sha256=part['probes_sha256'])
        total+=len(census);eligible+=sum(multiplicity.values());counts.update(c['source_disposition'] for c in census);all_probes.extend(probes)
    assert total==source['qualification_completion']['unique_designs']*2*len(fit['methods'])*len(fit['loading_modes'])*len(fit['trees'])
    planning=estimate(all_probes,fit);assert planning==json.loads((root/'conditional_planning.json').read_text())
    summary=dict(logical_cases=len(source['ids']),model_setting_rows=source['qualification_completion']['model_setting_rows'],
        unique_cohorts=len(manifest),candidate_rows=total,eligible_candidates=eligible,timing_groups=len(all_probes),
        source_status_counts=dict(counts),timing_status_counts=dict(Counter(p['status'] for p in all_probes)),
        measured_candidate_coverage=planning['measured_candidate_coverage'],unmeasured_review_candidate_coverage=planning['unmeasured_review_candidate_coverage'],
        conditional_budget_weighted_seconds=planning['conditional_budget_weighted_seconds'])
    assert all(receipt[k]==v for k,v in summary.items()) and receipt['production_finish_eta'] is None
    verify(bindings)
    result=dict(status='passed_full_shared_entity_timing_census_and_planning_accounting_readback',plan_sha256=sha(path),
        producer_receipt_sha256=sha(rp),fit_contract=source['fit_contract'],**summary,source_hashes=bindings,
        independent_numeric_or_hardware_timing_reimplementation=False,production_finish_eta=None,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.plan,a.output)
