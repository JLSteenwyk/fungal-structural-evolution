"""I/O and schemas for the complete target/background directed measurement catalog."""
import csv
import json
from collections import Counter
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
    # The closed usable-order summary binds the original target receipts, even
    # though downstream matching retained only its summary-level proof map.
    # Expand the full original receipt -> manifest -> checkpoint lineage here.
    for root,key in [(native,'native_receipt_sha256'),(diagnostic,'diagnostic_receipt_sha256'),(geometry,'geometry_receipt_sha256')]:bind(bindings,root/'receipt.json',sr[key])
    nr,dr,gr=[json.loads((root/'receipt.json').read_text()) for root in [native,diagnostic,geometry]]
    assert nr['status']=='complete_duplication_alignment_dispositions_pending_readback' and nr['plan_sha256']==sha(plan['target_native_plan'])
    assert nr['distinct_model_pairs']==sr['pairs'] and nr['directed_dispositions']==4*sr['pairs']
    assert dr['status']=='complete_primary_alignment_rmsd_diagnostic_not_scientific_acceptance' and dr['producer_receipt_sha256']==sha(native/'receipt.json')
    assert dr['counts']==nr['counts'] and dr['directed_dispositions']==nr['directed_dispositions']
    assert gr['status']=='complete_primary_diagnostic_geometry_pending_independent_readback' and gr['diagnostic_receipt_sha256']==sha(diagnostic/'receipt.json')
    assert gr['alignments']==dr['numerically_checked_alignments']==plan['expected']['target_numeric_states']
    for root,record in [(native,nr),(diagnostic,dr),(geometry,gr)]:
        for name,digest in record['artifacts'].items():bind(bindings,root/name,digest)
    native_keys=set();native_counts=Counter()
    with paths['target_checkpoints'].open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            pair,mask,os=Path(row['path']).stem.rsplit('-',2);order=int(os);key=pair,mask,order
            assert key not in native_keys and mask in ['full','plddt70'] and order in [0,1];native_keys.add(key);native_counts[mask+':'+row['status']]+=1
            assert row['path']==f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json'
            bind(bindings,native/row['path'],row['sha256'])
    native_pairs={p for p,m,o in native_keys}
    assert len(native_pairs)==sr['pairs'] and native_keys=={(p,m,o) for p in native_pairs for m in ['full','plddt70'] for o in [0,1]}
    assert dict(native_counts)==nr['counts']
    assert plan['screens']==mp['screens']==tp['screens']==bp['screens'] and c['selected_records']==plan['expected']['selected_records']
    for p in [plan['matching_plan'],plan['target_native_plan'],plan['target_summary_plan'],plan['target_summary_readback'],summary/'receipt.json',bp['measurement_completion'],union/'receipt.json',*paths.values()]:
        assert str(p) in bindings,('Source absent from full closed matching proof',str(p))
    verify(bindings)
    return dict(**paths,target_native_root=native,target_summary_receipt=sr,target_summary_readback=sa),bindings
