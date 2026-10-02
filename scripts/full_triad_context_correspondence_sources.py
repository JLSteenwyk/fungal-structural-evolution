"""Closed full correspondence comparison and original logical context sources."""
import json
from pathlib import Path
from full_triad_context_geometry_sources import GATES,POLICIES
from full_triad_correspondence_comparison_sources import METHODS,PERMUTATIONS,DEFINITIONS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

MASKS=['full','plddt70','both_masks']
METHOD_SCENARIOS=METHODS+['both_methods']
CORE_SCENARIOS=DEFINITIONS+['both_cores']
FULL_BITS=(1<<48)-1
SUMMARY_FIELDS=['target_contexts','context_design_records','reference_tie_records','duplicate_reference_links',
    'logical_reference_screen_decisions','context_screen_states','context_policy_decisions','summary_rows',
    'guide_contexts','reference_measurement_presence_counts','measured_triads','measured_comparison_groups']


def load_sources(plan,plan_path):
    bindings=dict(plan['pins']);bind(bindings,plan_path)
    cp=plan['comparison_completion'];c=json.loads(Path(cp).read_text())
    assert c['status']=='complete_verified_full_triad_correspondence_comparison' and c['exact_process_journals_checked']==2
    bind(bindings,cp);bind(bindings,c['full_hash_archive'],c['full_hash_archive_sha256'])
    archive=json.loads(Path(c['full_hash_archive']).read_text())
    assert archive['status']==c['status']+'_archive' and len(archive['services'])==2
    assert len(archive['source_hashes'])==c['bound_source_hashes']
    for q,d in archive['source_hashes'].items():bind(bindings,q,d)
    rp=Path(c['producer_receipt']);ap=Path(c['independent_readback'])
    bind(bindings,rp,c['producer_receipt_sha256']);bind(bindings,ap,c['independent_readback_sha256'])
    r,a=[json.loads(p.read_text()) for p in [rp,ap]]
    assert r['status']=='complete_full_triad_correspondence_comparison_pending_independent_readback'
    assert a['status']=='passed_full_triad_correspondence_raw_msa_sql_readback' and a['producer_receipt_sha256']==sha(rp)
    assert r['plan_sha256']==a['plan_sha256']==sha(plan['comparison_plan'])
    assert r['scientific_eligibility'] is a['scientific_eligibility'] is False
    for k,v in archive['summary'].items():assert r[k]==a[k]==c[k]==v
    assert c['source_ready_triads']==plan['expected']['measured_triads']==27056
    assert c['comparison_groups']==plan['expected']['measured_comparison_groups']==216448 and c['pair_states']==10389504
    assert plan['expected']==dict(target_contexts=283409,reference_tie_records=214461,duplicate_reference_links=428922,
        measured_triads=27056,measured_comparison_groups=216448,logical_reference_screen_decisions=34742682,
        context_screen_states=91824516,context_policy_decisions=459122580,summary_rows=3240)
    config=json.loads(Path(plan['comparison_plan']).read_text());assert config['screens']==plan['screens']
    for name,digest in r['artifacts'].items():bind(bindings,rp.parent/name,digest)
    wc=json.loads(Path(plan['triad_work_completion']).read_text())
    assert wc['status']=='complete_verified_full_reference_triad_work_design' and len(wc['services'])==2
    for path,digest in wc['source_hashes'].items():bind(bindings,path,digest)
    for key in ['target_contexts','reference_tie_records','duplicate_reference_links']:assert wc['summary'][key]==plan['expected'][key]
    for path in [plan['comparison_plan'],plan['triad_work_plan'],plan['triad_work_completion']]:bind(bindings,path)
    work=json.loads(Path(plan['triad_work_plan']).read_text());contexts=Path(work['output'])/'context_triad_design.jsonl.gz'
    assert str(contexts) in archive['source_hashes']
    verify(bindings)
    return contexts,rp.parent/'correspondence_comparison_groups.jsonl.gz',wc['summary'],bindings


def leaf_selectors(mask,method,core):
    return (MASKS[:2] if mask=='both_masks' else [mask],
            METHODS if method=='both_methods' else [method],
            DEFINITIONS if core=='both_cores' else [core])
