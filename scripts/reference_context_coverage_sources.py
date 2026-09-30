"""Closed full reference context and order/coverage source proof I/O."""
import json
from pathlib import Path
from reference_measurement_union_sources import bind, verify, artifacts
from run_ortholog_pair_guide_comparison import sha

POLICIES = ['lexical_pair_coverage', 'lexical_pair_coverage_native_own',
            'lexical_pair_coverage_native_both_guides', 'any_tie_pair_coverage_native_both_guides',
            'all_ties_pair_coverage_native_both_guides']


def load_sources(plan):
    bindings = dict(plan['pins'])
    config = json.loads(Path(plan['context_plan']).read_text()); root = Path(config['output'])
    r = json.loads((root / 'receipt.json').read_text()); a = json.loads(Path(plan['context_readback']).read_text())
    c = json.loads(Path(plan['context_completion']).read_text())
    assert r['status'] == 'complete_full_reference_context_measurement_design_pending_independent_readback'
    assert a['status'] == 'passed_full_reference_context_measurement_design_sql_readback'
    assert c['status'] == 'complete_verified_full_reference_context_measurement_design' and len(c['services']) == 2
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['context_plan'])
    assert a['producer_receipt_sha256'] == sha(root / 'receipt.json') == c['source_hashes'][str(root / 'receipt.json')]
    assert c['source_hashes'][plan['context_readback']] == sha(plan['context_readback'])
    assert all(r[k] == a[k] == v for k, v in c['summary'].items())
    for field in ['target_contexts', 'reference_tie_records', 'duplicate_reference_links']:
        assert r[field] == plan['expected'][field]
    assert r['availability_side_links_checked'] == plan['expected']['availability_side_links']
    artifacts(bindings, root, r)
    assert c['source_hashes'][str(root / 'context_measurement_design.jsonl')] == r['artifacts']['context_measurement_design.jsonl']
    coverage_config = json.loads(Path(plan['coverage_plan']).read_text()); coverage_root = Path(coverage_config['output'])
    cr = json.loads((coverage_root / 'receipt.json').read_text()); cc = json.loads(Path(plan['coverage_completion']).read_text())
    assert cc['status'] == 'complete_verified_full_reference_order_and_original_coverage' and cc['exact_process_journals_checked'] == 2
    ca = json.loads(Path(cc['independent_readback']).read_text())
    assert cr['status'] == 'complete_full_reference_order_and_original_coverage_pending_independent_readback'
    assert ca['status'] == 'passed_full_reference_order_and_original_coverage_readback'
    assert cr['plan_sha256'] == ca['plan_sha256'] == sha(plan['coverage_plan'])
    assert cc['producer_receipt'] == str(coverage_root / 'receipt.json')
    assert cc['producer_receipt_sha256'] == ca['producer_receipt_sha256'] == sha(coverage_root / 'receipt.json')
    bind(bindings, cc['independent_readback'], cc['independent_readback_sha256'])
    bind(bindings, cc['full_hash_archive'], cc['full_hash_archive_sha256'])
    archive = json.loads(Path(cc['full_hash_archive']).read_text()); assert len(archive['services']) == 2
    fields = ['full_pairs', 'directed_dispositions', 'pair_mask_rows', 'pair_screen_decisions', 'order_summary_counts',
              'maximum_order_differences', 'pair_pass_counts', 'pair_exclusion_counts', 'pair_both_masks_pass_counts',
              'source_dispositions', 'native_status_counts', 'numerical_counts']
    for field in fields: assert cr[field] == ca[field] == cc[field] == archive['summary'][field]
    assert cr['full_pairs'] == plan['expected']['model_pairs'] and cr['pair_mask_rows'] == 2 * cr['full_pairs']
    assert cr['screens'] == coverage_config['screens'] == plan['screens'] and len(plan['screens']) == 6
    artifacts(bindings, coverage_root, cr)
    assert archive['source_hashes'][str(coverage_root / 'pair_mask_order_coverage.tsv')] == cr['artifacts']['pair_mask_order_coverage.tsv']
    union_config = json.loads(Path(coverage_config['union_plan']).read_text())
    native_config = json.loads(Path(union_config['native_plan']).read_text())
    assert native_config['inventory'] == config['inventory'] and native_config['base_queue'] == config['queue']
    for path in [plan['context_plan'], plan['context_readback'], plan['context_completion'], plan['coverage_plan'],
                 plan['coverage_completion'], coverage_config['union_plan'], union_config['native_plan']]: bind(bindings, path)
    verify(bindings)
    return root / 'context_measurement_design.jsonl', coverage_root / 'pair_mask_order_coverage.tsv', r, cr, bindings
