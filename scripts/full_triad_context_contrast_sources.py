"""Closed contrast sensitivities and unchanged original context source lineage."""
import json
from pathlib import Path
from full_triad_contrast_sensitivity_sources import MASKS,CORES
from full_triad_context_geometry_sources import GATES,POLICIES
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

PHYSICAL_DIRECTIONS=['positive','negative','within_numerical_tolerance','sign_uncertain']
STATES=PHYSICAL_DIRECTIONS+['reference_direction_disagreement','no_tied_reference','source_gate_excluded',
                          'geometry_quality_excluded','no_eligible_reference','not_all_ties_eligible']
SUMMARY_FIELDS=['target_contexts','context_design_records','reference_tie_records','duplicate_reference_links',
    'logical_reference_screen_decisions','context_screen_states','context_policy_decisions','summary_rows',
    'guide_contexts','reference_measurement_presence_counts','measured_triads','measured_contrast_groups','direction_counts']


def load_sources(plan,path):
    bindings=dict(plan['pins']);bind(bindings,path)
    cp=plan['contrast_completion'];c=json.loads(Path(cp).read_text());bind(bindings,cp)
    assert c['status']=='complete_verified_full_triad_contrast_sensitivity' and c['exact_process_journals_checked']==2
    bind(bindings,c['full_hash_archive'],c['full_hash_archive_sha256']);archive=json.loads(Path(c['full_hash_archive']).read_text())
    assert archive['status']==c['status']+'_archive' and len(archive['services'])==2 and len(archive['source_hashes'])==c['bound_source_hashes']
    for q,d in archive['source_hashes'].items():bind(bindings,q,d)
    rp=Path(c['producer_receipt']);ap=Path(c['independent_readback'])
    bind(bindings,rp,c['producer_receipt_sha256']);bind(bindings,ap,c['independent_readback_sha256'])
    r,a=[json.loads(p.read_text()) for p in [rp,ap]]
    assert r['status']=='complete_full_triad_contrast_sensitivity_pending_independent_readback'
    assert a['status']=='passed_full_triad_contrast_sensitivity_raw_fit_sql_readback' and a['producer_receipt_sha256']==sha(rp)
    assert r['plan_sha256']==a['plan_sha256']==sha(plan['contrast_plan']) and r['scientific_eligibility'] is a['scientific_eligibility'] is False
    for k,v in archive['summary'].items():assert c[k]==r[k]==a[k]==v
    assert c['measured_triads']==plan['expected']['measured_triads']==27056 and c['sensitivity_groups']==plan['expected']['measured_contrast_groups']==243504
    config=json.loads(Path(plan['contrast_plan']).read_text());assert config['screens']==plan['screens'] and config['contrast_numerical_tolerance_angstrom']==1e-9
    assert config['triad_work_plan']==plan['triad_work_plan'] and config['triad_work_completion']==plan['triad_work_completion']
    for name,d in r['artifacts'].items():bind(bindings,rp.parent/name,d)
    wc=json.loads(Path(plan['triad_work_completion']).read_text())
    assert wc['status']=='complete_verified_full_reference_triad_work_design' and len(wc['services'])==2
    for q,d in wc['source_hashes'].items():bind(bindings,q,d)
    for k in ['target_contexts','reference_tie_records','duplicate_reference_links']:assert wc['summary'][k]==plan['expected'][k]
    assert plan['expected']['logical_reference_screen_decisions']==11580894 and plan['expected']['context_screen_states']==30608172
    assert plan['expected']['context_policy_decisions']==153040860 and plan['expected']['summary_rows']==10800
    work=json.loads(Path(plan['triad_work_plan']).read_text());contexts=Path(work['output'])/'context_triad_design.jsonl.gz'
    assert str(contexts) in archive['source_hashes']
    for p in [plan['contrast_plan'],plan['triad_work_plan'],plan['triad_work_completion']]:bind(bindings,p)
    verify(bindings)
    return contexts,rp.parent/'contrast_sensitivity.jsonl.gz',wc['summary'],bindings
