"""Closed case-index/catalog sources and schemas for all matched measurements."""
import json
from pathlib import Path
from background_measurement_union_sources import closed_source
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

OUTCOMES=['rmsd','native_tm_dissimilarity']
STATE_FIELDS=['state_key','native_status','numerical_usable','numerical_exclusion_reasons','aligned_length','rmsd','native_tm_dissimilarity']
STATES=[(r,o) for r in ['target','background'] for o in [0,1]]
BASE_FIELDS=['row_identity','case_id','physical_case_id','target_id','background_id','guide','mask','target_sequence_distance','background_sequence_distance','target_same_model','background_same_model','target_pair_key','background_pair_key','target_mask_pass_bits','control_mask_pass_bits','joint_mask_pass_bits','target_both_pass_bits','control_both_pass_bits','joint_both_pass_bits']
FIELDS=BASE_FIELDS+[f'{r}_order{o}_{f}' for r,o in STATES for f in STATE_FIELDS]+[f'{outcome}_{f}' for outcome in OUTCOMES for f in ['target_both_orders_mean','background_both_orders_mean','both_orders_mean_delta','order_pair_00_delta','order_pair_01_delta','order_pair_10_delta','order_pair_11_delta','usable_order_pair_bits','complete_order_envelope_min','complete_order_envelope_max','complete_order_envelope_span','disposition']]
SUMMARY_FIELDS=['logical_cases','physical_cases','selected_records','unmatched_decisions','case_mask_rows','order_pair_cells','valid_order_pair_cells','complete_envelopes','incomplete_envelopes','same_model_case_mask_rows','outcomes','masks']


def load(plan,path):
    bindings=dict(plan['pins']);bind(bindings,path);sources={}
    for key,status in [('case_index','complete_verified_full_matching_logical_case_index'),('catalog','complete_verified_full_expanded_directed_measurement_catalog')]:
        c=closed_source(plan[key+'_completion'],status,status+'_archive',2,bindings);assert c['scientific_eligibility'] is False
        config=json.loads(Path(plan[key+'_plan']).read_text());root=Path(config['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
        assert c['producer_receipt']==str(rp) and bindings[str(rp)]==c['producer_receipt_sha256'] and r['plan_sha256']==sha(plan[key+'_plan'])
        assert r['scientific_eligibility'] is False
        for n,d in r['artifacts'].items():bind(bindings,root/n,d)
        sources[key]=dict(root=root,receipt=r,completion=c,config=config)
    case,catalog=sources['case_index'],sources['catalog']
    assert case['config']['matching_completion']==catalog['config']['matching_completion']==plan['matching_completion']
    assert case['config']['matching_plan']==catalog['config']['matching_plan']
    assert case['receipt']['selected_records']==plan['expected']['selected_records']
    assert case['receipt']['logical_cases']==plan['expected']['logical_cases']
    assert case['receipt']['unmatched_decisions']==plan['expected']['unmatched_decisions']
    assert catalog['receipt']['directed_states']==plan['expected']['directed_states']
    assert catalog['receipt']['masks']==['full','plddt70']
    verify(bindings);return sources,bindings
