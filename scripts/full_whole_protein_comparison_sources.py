"""Require exhaustive closed optimization evidence before full comparison integration."""
import json
from pathlib import Path
from whole_protein_support_registry import SupportRegistry
from close_full_whole_protein_optimization import index
from reference_measurement_union_sources import bind,verify
from screen_duplication_alignment_reuse import sha

SUMMARY_FIELDS=['fit_summaries','parameter_records','comparisons','trees','original_fit_status_counts',
    'effective_fit_status_counts','followup_status_counts','selected_source_kind_counts',
    'comparison_status_counts','relation_counts','followup_cases_integrated']


def load_sources(plan,path):
    bindings=dict(plan['pins']);bind(bindings,path)
    cp=Path(plan['optimization_completion']);c=json.loads(cp.read_text());bind(bindings,cp)
    expected_journals=plan['expected_process_journals']
    historical=plan.get('historical_source_hashes',{})
    assert type(expected_journals) is int and expected_journals>0
    assert c['status']=='complete_verified_full_whole_protein_optimization' and c['scientific_eligibility'] is False and c['exact_process_journals_checked']==expected_journals
    assert c['historical_source_hashes']==len(historical)
    bind(bindings,c['full_hash_archive'],c['full_hash_archive_sha256']);a=json.loads(Path(c['full_hash_archive']).read_text())
    assert a['status']==c['status']+'_archive' and len(a['services'])==expected_journals and len(a['source_hashes'])==c['bound_source_hashes']
    assert a.get('historical_source_hashes',{})==historical
    for k,v in a['summary'].items():assert c[k]==v
    for p,h in a['source_hashes'].items():
        assert p not in historical
        bind(bindings,p,h)
    assert c['unique_inputs']==75070 and c['full_dispositions']==375350
    op=json.loads(Path(plan['optimization_plan']).read_text());assert c['completion_plan_sha256']==sha(plan['optimization_plan'])
    assert op['fit_plan']==plan['fit_plan'] and op['followup_plan']==plan['followup_plan']
    fp=json.loads(Path(plan['fit_plan']).read_text());pp=json.loads(Path(plan['followup_plan']).read_text())
    roots=[Path(fp['output']),Path(pp['output'])];original=index(roots[0]/'fit_manifest.jsonl');follow=index(roots[1]/'followup_manifest.jsonl')
    assert len(original)==375350 and len(follow)==c['flagged_fits']
    assert set(follow)=={k for k,r in original.items() if r['status']=='ml_candidate_requires_optimization_review'}
    assert len({k[0] for k in original})==75070 and len({k[1] for k in original})==5
    for row in original.values():assert bindings[row['path']]==row['sha256']
    for row in follow.values():assert bindings[str(roots[1]/row['path'])]==row['sha256']
    links=Path(plan['links']);lp=json.loads(Path(plan['link_proof']).read_text());lr=json.loads((links/'receipt.json').read_text())
    assert lp['status']=='passed_full_whole_protein_comparison_input_link_readback' and lp['source_receipt_sha256']==sha(links/'receipt.json')
    for p,h in lp.get('source_hashes',{}).items():bind(bindings,p,h)
    bind(bindings,plan['link_proof']);bind(bindings,links/'receipt.json')
    for n,h in lr['artifacts'].items():bind(bindings,links/n,h)
    support=SupportRegistry(plan['support'],plan['support_proof']);bind(bindings,Path(plan['support'])/'receipt.json',support.receipt_sha256);bind(bindings,plan['support_proof'],support.proof_sha256)
    sp=json.loads(Path(plan['support_proof']).read_text())
    for p,h in sp['source_hashes'].items():bind(bindings,p,h)
    for n,h in json.loads((Path(plan['support'])/'receipt.json').read_text())['artifacts'].items():bind(bindings,Path(plan['support'])/n,h)
    assert plan['expected']==dict(unique_inputs=75070,fit_summaries=375350,parameter_records=375350,trees=5,comparisons_per_tree=829440,comparisons=4147200)
    verify(bindings);return original,follow,roots[1],lr,support,bindings
