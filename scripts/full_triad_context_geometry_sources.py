"""Closed all-order results and full unchanged context work-design source I/O."""
import json
from pathlib import Path
from full_triad_robustness_sources import DEFINITIONS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

MASKS = ['full', 'plddt70', 'both_masks']
GATES = ['source_design_ready_for_correspondence', 'source_design_ready_native_own', 'source_design_ready_native_both_guides']
POLICIES = ['lexical_geometry', 'lexical_geometry_native_own', 'lexical_geometry_native_both_guides',
            'any_tie_geometry_native_both_guides', 'all_ties_geometry_native_both_guides']
SUMMARY_FIELDS = ['target_contexts', 'context_design_records', 'reference_tie_records', 'duplicate_reference_links',
                  'logical_reference_screen_decisions', 'context_screen_states', 'context_policy_decisions', 'summary_rows',
                  'guide_contexts', 'reference_measurement_presence_counts', 'measured_triads', 'measured_robustness_groups']


def load_sources(plan, plan_path):
    bindings = dict(plan['pins']); bind(bindings, plan_path)
    c = json.loads(Path(plan['robustness_completion']).read_text()); bind(bindings, plan['robustness_completion'])
    assert c['status'] == 'complete_verified_full_triad_all_order_robustness' and c['exact_process_journals_checked'] == 2
    bind(bindings, c['full_hash_archive'], c['full_hash_archive_sha256'])
    archive = json.loads(Path(c['full_hash_archive']).read_text())
    assert archive['status'] == 'complete_verified_full_triad_all_order_robustness_archive' and len(archive['services']) == 2
    assert len(archive['source_hashes']) == c['bound_source_hashes']
    for path, digest in archive['source_hashes'].items(): bind(bindings, path, digest)
    rp, ap = Path(c['producer_receipt']), Path(c['independent_readback'])
    bind(bindings, rp, c['producer_receipt_sha256']); bind(bindings, ap, c['independent_readback_sha256'])
    r, a = [json.loads(p.read_text()) for p in [rp, ap]]
    assert r['status'] == 'complete_full_triad_all_order_robustness_pending_independent_readback'
    assert a['status'] == 'passed_full_triad_all_order_robustness_sql_readback' and a['producer_receipt_sha256'] == sha(rp)
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['robustness_plan'])
    assert r['correspondence_work_triads'] == a['correspondence_work_triads'] == c['correspondence_work_triads'] == plan['expected']['measured_triads']
    assert r['robustness_groups'] == a['robustness_groups'] == c['robustness_groups'] == plan['expected']['measured_robustness_groups']
    config = json.loads(Path(plan['robustness_plan']).read_text()); assert config['screens'] == plan['screens']
    assert config['triad_work_plan'] == plan['triad_work_plan'] and config['triad_work_completion'] == plan['triad_work_completion']
    for name, digest in r['artifacts'].items(): bind(bindings, rp.parent / name, digest)
    wc = json.loads(Path(plan['triad_work_completion']).read_text())
    assert wc['status'] == 'complete_verified_full_reference_triad_work_design' and len(wc['services']) == 2
    for path, digest in wc['source_hashes'].items(): bind(bindings, path, digest)
    for key in ['target_contexts', 'reference_tie_records', 'duplicate_reference_links']: assert wc['summary'][key] == plan['expected'][key]
    for path in [plan['robustness_plan'], plan['triad_work_plan'], plan['triad_work_completion']]: bind(bindings, path)
    work = json.loads(Path(plan['triad_work_plan']).read_text()); contexts = Path(work['output']) / 'context_triad_design.jsonl.gz'
    assert str(contexts) in bindings
    verify(bindings)
    return contexts, rp.parent / 'triad_order_robustness.jsonl.gz', wc['summary'], bindings
