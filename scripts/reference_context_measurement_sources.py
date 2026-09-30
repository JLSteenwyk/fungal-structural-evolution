"""Closed native-reference context and frozen gene/model ledger source I/O."""
import csv
import hashlib
import json
from pathlib import Path
from reference_measurement_union_sources import bind, verify, artifacts, pair_ends
from run_ortholog_pair_guide_comparison import sha


def event_key(guide, row):
    return (guide, row['family'], row['gene_node'], *sorted([row['gene_a'], row['gene_b']]))


def pair_hash(ends):
    return hashlib.sha256(json.dumps(sorted(ends), separators=(',', ':')).encode()).hexdigest()


def load_sources(plan):
    bindings = dict(plan['pins']); root = Path(plan['native_contexts'])
    rp = root / 'receipt.json'; r = json.loads(rp.read_text())
    proof = json.loads(Path(plan['native_readback']).read_text()); closed = json.loads(Path(plan['native_completion']).read_text())
    assert r['status'] == 'complete_full_reference_native_orthology_pending_readback'
    assert proof['status'] == 'passed_full_reference_native_orthology_readback'
    assert closed['status'] == 'complete_verified_full_reference_native_orthology' and len(closed['services']) == 2
    assert proof['producer_receipt_sha256'] == sha(rp) == closed['source_hashes'][str(rp)]
    assert closed['source_hashes'][plan['native_readback']] == sha(plan['native_readback'])
    for field in ['target_contexts', 'reference_tie_records', 'duplicate_reference_links']:
        assert r[field] == proof[field] == closed['summary'][field] == plan['expected'][field]
    assert r['guide_contexts'] == closed['summary']['guide_contexts']
    artifacts(bindings, root, r)
    for name in r['artifacts']:
        assert closed['source_hashes'][str(root / name)] == r['artifacts'][name]
    for path in [plan['native_readback'], plan['native_completion']]: bind(bindings, path)
    queue = Path(plan['queue']); qr = json.loads((queue / 'receipt.json').read_text())
    qp = json.loads(Path(plan['queue_readback']).read_text())
    assert qr['status'] == 'complete_reviewed_duplication_model_pair_queue'
    assert qp['status'] == 'passed_full_duplication_model_pair_queue_export_readback'
    assert qr['counts']['event_links'] == plan['expected']['target_contexts']
    for name in ['receipt.json', 'event_model_pair_links.tsv', 'models.jsonl', 'model_pairs.tsv']:
        p = queue / name; assert qp['source_sha256'][str(p)] == sha(p)
        if name != 'receipt.json': assert qr['artifacts'][name] == qp['source_sha256'][str(p)]
        bind(bindings, p, qp['source_sha256'][str(p)])
    bind(bindings, plan['queue_readback'])
    inventory = Path(plan['inventory']); ir = json.loads((inventory / 'receipt.json').read_text())
    ip = json.loads(Path(plan['inventory_readback']).read_text())
    ic = json.loads(Path(plan['inventory_completion']).read_text())
    assert ir['status'] == 'complete_provisional_reference_comparison_inventory'
    assert ip['status'] == 'passed_full_reference_comparison_ledger_and_native_model_readback'
    assert ic['status'] == 'complete_verified_expanded_duplication_reference_pipeline' and len(ic['services']) == 6
    assert ip['producer_receipt_sha256'] == sha(inventory / 'receipt.json')
    assert ic['source_hashes'][str(inventory / 'receipt.json')] == ip['producer_receipt_sha256']
    assert ic['source_hashes'][plan['inventory_readback']] == sha(plan['inventory_readback'])
    assert ip['event_reference_comparisons'] == ic['summary']['event_reference_comparisons'] == plan['expected']['availability_side_links']
    assert ir['unique_distinct_model_pairs'] == ip['distinct_model_pairs'] == plan['expected']['model_pairs']
    artifacts(bindings, inventory, ir)
    for path in [plan['inventory_readback'], plan['inventory_completion']]: bind(bindings, path)
    # The full independent reference ledger already binds the actual primary queue.
    assert ip['primary_queue_readback_sha256'] == sha(plan['queue_readback'])
    assert ir['source_pins'][str(queue / 'receipt.json')] == sha(queue / 'receipt.json')
    verify(bindings)
    return root / 'contexts.jsonl', queue, inventory, r, bindings


def load_tables(queue, inventory):
    duplicates, pairs, links = {}, {}, {}
    with (queue / 'event_model_pair_links.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = event_key(row['guide'], row); assert key not in duplicates
            duplicates[key] = row
    with (inventory / 'model_pairs.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            ends = pair_ends(row); assert ends[0] != ends[1] and row['pair_key'] == pair_hash(ends)
            assert row['pair_key'] not in pairs; pairs[row['pair_key']] = row
    with (inventory / 'event_reference_comparisons.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = (*event_key(row['guide'], row), row['reference_gene'], row['focal_gene'])
            assert key not in links and row['focal_gene'] == row['gene_' + row['focal_side']]
            links[key] = row
    return duplicates, pairs, links
