"""I/O and schemas for the complete target/background directed measurement catalog."""
import json
from pathlib import Path
from background_measurement_union_sources import closed_source,NUMERIC_INTS,NUMERIC_FLOATS,GEOMETRY_INTS,GEOMETRY_FLOATS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

NUMERIC_FIELDS=NUMERIC_INTS+NUMERIC_FLOATS+['rmsd_status']
GEOMETRY_FIELDS=[n for n in GEOMETRY_INTS+GEOMETRY_FLOATS if n!='aligned_length']+['geometry_status']
FIELDS=['role','pair_key','mask','order','model_a','version_a','model_b','version_b','original_length_a','original_length_b',
    'source_kind','source_checkpoint','source_checkpoint_sha256','source_order','native_status','numerical_usable','numerical_exclusion_reasons',
    'pair_mask_pass_bits','both_masks_pass_bits','original_coverage_a','original_coverage_b']+NUMERIC_FIELDS+GEOMETRY_FIELDS
SUMMARY_FIELDS=['target_pairs','background_pairs','pair_mask_rows','directed_states','numeric_states','usable_states','unavailable_states','native_status_counts','numerical_exclusion_counts','source_kind_counts','screens','masks']


def load(plan,path):
    bindings=dict(plan['pins']);bind(bindings,path)
    c=closed_source(plan['matching_completion'],'complete_verified_full_fixed_matched_coverage_attrition','complete_verified_full_matched_coverage_archive',2,bindings)
    assert c['scientific_eligibility'] is False
    mp=json.loads(Path(plan['matching_plan']).read_text());assert c['source_plan']==plan['matching_plan'] and c['source_plan_sha256']==sha(plan['matching_plan'])
    tp=json.loads(Path(mp['target_plan']).read_text());bp=json.loads(Path(mp['background_plan']).read_text())
    summary=Path(tp['summary']);sp=json.loads(Path(plan['target_summary_plan']).read_text())
    assert Path(sp['output'])==summary and tp['summary_audit']==plan['target_summary_readback']
    sr,sa=[json.loads(p.read_text()) for p in [summary/'receipt.json',Path(plan['target_summary_readback'])]]
    assert sr['status']=='complete_primary_numerically_usable_order_summary_pending_independent_readback' and sa['status']=='passed_full_primary_usable_order_summary_readback'
    assert sr['plan_sha256']==sha(plan['target_summary_plan']) and sa['source_receipt_sha256']==sha(summary/'receipt.json')
    assert sr['pairs']==plan['expected']['target_pairs'] and sa['pair_mask_rows']==sr['pair_mask_rows']==2*sr['pairs']
    up=json.loads(Path(bp['measurement_plan']).read_text());union=Path(up['output'])
    uc=json.loads(Path(bp['measurement_completion']).read_text());assert uc['status']=='complete_verified_full_background_measurement_union' and uc['full_pairs']==plan['expected']['background_pairs']
    assert uc['producer_receipt']==str(union/'receipt.json') and sha(union/'receipt.json')==uc['producer_receipt_sha256']
    ur=json.loads((union/'receipt.json').read_text());assert ur['directed_dispositions']==4*uc['full_pairs'] and ur['status']=='complete_full_background_measurement_union_pending_independent_readback'
    native,diagnostic,geometry=[Path(sp[k]) for k in ['native','diagnostic','geometry']]
    np=json.loads(Path(plan['target_native_plan']).read_text());assert Path(np['output'])==native
    paths={'target_quality':Path(tp['output'])/'pair_mask_coverage.tsv','background_quality':Path(bp['output'])/'pair_mask_coverage.tsv','target_summary':summary/'pair_mask_order_summary.tsv','target_checkpoints':native/'checkpoint_manifest.tsv','target_pairs':Path(np['queue'])/'model_pairs.tsv','target_numeric':diagnostic/'numeric_readback.tsv','target_geometry':geometry/'alignment_geometry.tsv','background_states':union/'background_measurement_dispositions.jsonl.gz'}
    assert plan['screens']==mp['screens']==tp['screens']==bp['screens'] and c['selected_records']==plan['expected']['selected_records']
    for p in [plan['matching_plan'],plan['target_native_plan'],plan['target_summary_plan'],plan['target_summary_readback'],summary/'receipt.json',bp['measurement_completion'],union/'receipt.json',*paths.values()]:
        assert str(p) in bindings,('Source absent from full closed matching proof',str(p))
    verify(bindings)
    return dict(**paths,target_native_root=native,target_summary_receipt=sr,target_summary_readback=sa),bindings
