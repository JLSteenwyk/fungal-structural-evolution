"""Closed full reference-context work design and original pair-ledger source I/O."""
import json
from pathlib import Path
from reference_measurement_union_sources import bind, verify, artifacts
from run_ortholog_pair_guide_comparison import sha

SUMMARY_FIELDS = ['target_contexts', 'context_design_records', 'reference_tie_records',
                  'duplicate_reference_links', 'availability_side_links_preserved', 'guide_contexts',
                  'unique_ordered_model_triads', 'correspondence_work_triads', 'potential_correspondence_states',
                  'triad_model_identity_counts', 'duplicate_design_counts', 'reference_link_counts', 'empty_reference_context_counts']


def load_sources(plan):
    bindings = dict(plan['pins']); config = json.loads(Path(plan['context_plan']).read_text())
    root = Path(config['output']); rp = root / 'receipt.json'
    r = json.loads(rp.read_text()); a = json.loads(Path(plan['context_readback']).read_text())
    c = json.loads(Path(plan['context_completion']).read_text())
    assert r['status'] == 'complete_full_reference_context_measurement_design_pending_independent_readback'
    assert a['status'] == 'passed_full_reference_context_measurement_design_sql_readback'
    assert c['status'] == 'complete_verified_full_reference_context_measurement_design' and len(c['services']) == 2
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['context_plan'])
    assert a['producer_receipt_sha256'] == sha(rp) == c['source_hashes'][str(rp)]
    assert c['source_hashes'][plan['context_readback']] == sha(plan['context_readback'])
    assert all(r[k] == a[k] == v for k, v in c['summary'].items())
    for key in ['target_contexts', 'reference_tie_records', 'duplicate_reference_links']:
        assert r[key] == plan['expected'][key]
    assert r['availability_side_links_checked'] == plan['expected']['availability_side_links']
    for path, h in c['source_hashes'].items(): bind(bindings, path, h)
    artifacts(bindings, root, r)
    queue, inventory = Path(config['queue']), Path(config['inventory'])
    qr = json.loads((queue / 'receipt.json').read_text()); qp = json.loads(Path(config['queue_readback']).read_text())
    ir = json.loads((inventory / 'receipt.json').read_text())
    assert qr['status'] == 'complete_reviewed_duplication_model_pair_queue'
    assert qp['status'] == 'passed_full_duplication_model_pair_queue_export_readback'
    assert ir['unique_distinct_model_pairs'] == plan['expected']['reference_pairs']
    for name in ['receipt.json', 'event_model_pair_links.tsv', 'model_pairs.tsv', 'models.jsonl']:
        path = str(queue / name)
        assert bindings[path] == qp['source_sha256'][path]
        if name != 'receipt.json': assert bindings[path] == qr['artifacts'][name]
    assert bindings[str(inventory / 'model_pairs.tsv')] == ir['artifacts']['model_pairs.tsv']
    for path in [plan['context_plan'], plan['context_readback'], plan['context_completion']]: bind(bindings, path)
    verify(bindings)
    return root / 'context_measurement_design.jsonl', queue, inventory, r, bindings
