"""Closed sequence/structural correspondences for complete joint direction sensitivity."""
import json
from pathlib import Path
from full_triad_correspondence_comparison_sources import METHODS,PERMUTATIONS,DEFINITIONS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

MASKS=['full','plddt70','both_masks']
METHOD_SCENARIOS=METHODS+['both_methods']
CORE_SCENARIOS=DEFINITIONS+['both_cores']
DIRECTIONS=['positive','negative','within_numerical_tolerance','sign_uncertain','unavailable','nonunique_fit']
QUALIFIED=DIRECTIONS+['excluded_by_quality']
FULL_BITS=(1<<48)-1
SUMMARY_FIELDS=['measured_triads','sequence_fit_rows','structural_fit_rows','source_sequence_groups','source_structural_groups',
    'source_comparison_groups','joint_direction_groups','screen_decisions','summary_rows','joint_direction_counts',
    'qualified_direction_counts','source_category_relation_counts','maximum_joint_contrast_span']


def bundle(path,status,producer_status,reader_status,plan,bindings):
    c=json.loads(Path(path).read_text());assert c['status']==status and c['exact_process_journals_checked']==2
    bind(bindings,path);bind(bindings,c['full_hash_archive'],c['full_hash_archive_sha256']);a=json.loads(Path(c['full_hash_archive']).read_text())
    assert a['status']==status+'_archive' and len(a['services'])==2 and len(a['source_hashes'])==c['bound_source_hashes']
    for p,h in a['source_hashes'].items():bind(bindings,p,h)
    rp=Path(c['producer_receipt']);ap=Path(c['independent_readback']);bind(bindings,rp,c['producer_receipt_sha256']);bind(bindings,ap,c['independent_readback_sha256'])
    r,v=[json.loads(p.read_text()) for p in [rp,ap]]
    assert r['status']==producer_status and v['status']==reader_status and v['producer_receipt_sha256']==sha(rp)
    assert r['plan_sha256']==v['plan_sha256']==sha(plan) and r['scientific_eligibility'] is v['scientific_eligibility'] is False
    for k,value in a['summary'].items():assert c[k]==r[k]==v[k]==value
    for p,h in r['artifacts'].items():bind(bindings,rp.parent/p,h)
    bind(bindings,plan);return c,rp.parent


def load_sources(plan,path):
    bindings=dict(plan['pins']);bind(bindings,path)
    comparison,root=bundle(plan['comparison_completion'],'complete_verified_full_triad_correspondence_comparison',
        'complete_full_triad_correspondence_comparison_pending_independent_readback','passed_full_triad_correspondence_raw_msa_sql_readback',plan['comparison_plan'],bindings)
    structural,sroot=bundle(plan['structural_robustness_completion'],'complete_verified_full_triad_all_order_robustness',
        'complete_full_triad_all_order_robustness_pending_independent_readback','passed_full_triad_all_order_robustness_sql_readback',plan['structural_robustness_plan'],bindings)
    assert comparison['source_ready_triads']==structural['correspondence_work_triads']==27056 and comparison['ordered_model_triads']==31235
    assert comparison['sequence_fit_rows']==649344 and comparison['structural_fit_rows']==structural['fit_rows']==865792
    assert comparison['sequence_groups']==structural['robustness_groups']==108224 and comparison['comparison_groups']==216448 and comparison['pair_states']==10389504
    cp=json.loads(Path(plan['comparison_plan']).read_text());sp=json.loads(Path(plan['structural_robustness_plan']).read_text())
    assert cp['structural_geometry_plan']==sp['geometry_plan'] and cp['screens']==sp['screens']==plan['screens']
    work=Path(sp['triad_work_plan']);bind(bindings,work);catalog=Path(json.loads(work.read_text())['output'])/'ordered_model_triads.jsonl'
    assert str(catalog) in bindings
    with catalog.open() as f:triads=[r for line in f if (r:=json.loads(line))['source_design_ready_links']]
    assert len(triads)==27056
    sequence=json.loads(Path(cp['sequence_geometry_plan']).read_text());geometry=json.loads(Path(cp['structural_geometry_plan']).read_text())
    raw_sequence=Path(sequence['output'])/'sequence_common_residue_fits.tsv.gz';raw_structural=Path(geometry['output'])/'common_residue_fits.tsv.gz'
    paths=dict(sequence_groups=root/'sequence_order_robustness.jsonl.gz',structural_groups=sroot/'triad_order_robustness.jsonl.gz',
               comparison_groups=root/'correspondence_comparison_groups.jsonl.gz',raw_sequence=raw_sequence,raw_structural=raw_structural)
    assert all(str(p) in bindings for p in paths.values())
    assert plan['expected']==dict(measured_triads=27056,sequence_fit_rows=649344,structural_fit_rows=865792,source_sequence_groups=108224,
        source_structural_groups=108224,source_comparison_groups=216448,joint_direction_groups=730512,screen_decisions=4383072,summary_rows=1134)
    assert plan['contrast_numerical_tolerance_angstrom']==1e-9
    verify(bindings);return paths,triads,bindings


def selectors(mask,method,core):
    return (MASKS[:2] if mask=='both_masks' else [mask],METHODS if method=='both_methods' else [method],DEFINITIONS if core=='both_cores' else [core])
