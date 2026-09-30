"""Closed full-reference union and original-length input/provenance I/O."""
import csv
import hashlib
import json
from collections import Counter
from itertools import groupby
from pathlib import Path

from reference_measurement_union_sources import bind, verify, artifacts, pair_ends
from run_ortholog_pair_guide_comparison import sha


def load_sources(plan):
    bindings = dict(plan['pins'])
    union_plan = json.loads(Path(plan['union_plan']).read_text())
    native_plan = json.loads(Path(union_plan['native_plan']).read_text())
    root = Path(union_plan['output'])
    closure = json.loads(Path(plan['union_completion']).read_text())
    assert closure['status'] == 'complete_verified_full_reference_measurement_union'
    assert closure['full_pairs'] == plan['expected']['pairs']
    assert closure['directed_dispositions'] == 4 * plan['expected']['pairs']
    assert closure['exact_process_journals_checked'] == 2
    bind(bindings, closure['full_hash_archive'], closure['full_hash_archive_sha256'])
    archive = json.loads(Path(closure['full_hash_archive']).read_text())
    assert len(archive['services']) == 2
    assert archive['summary']['full_pairs'] == closure['full_pairs']
    assert archive['summary']['directed_dispositions'] == closure['directed_dispositions']
    assert closure['producer_receipt'] == str(root / 'receipt.json')
    receipt = json.loads((root / 'receipt.json').read_text())
    proof = json.loads(Path(closure['independent_readback']).read_text())
    assert receipt['status'] == 'complete_full_reference_measurement_union_pending_independent_readback'
    assert proof['status'] == 'passed_full_reference_measurement_union_readback'
    assert receipt['plan_sha256'] == proof['plan_sha256'] == sha(plan['union_plan'])
    assert closure['producer_receipt_sha256'] == proof['producer_receipt_sha256'] == sha(root / 'receipt.json')
    bind(bindings, closure['independent_readback'], closure['independent_readback_sha256'])
    for field in ['full_pairs', 'directed_dispositions', 'source_dispositions', 'native_status_counts', 'numerical_counts']:
        assert receipt[field] == proof[field] == closure[field] == archive['summary'][field]
    artifacts(bindings, root, receipt)
    assert archive['source_hashes'][str(root / 'reference_measurement_dispositions.jsonl')] == receipt['artifacts']['reference_measurement_dispositions.jsonl']
    pairs, models, inputs = load_inventory_inputs(native_plan, plan, bindings)
    screens = json.loads(Path(plan['screens_plan']).read_text())['screens']
    assert screens == plan['screens'] and len(screens) == 6 and len({s['id'] for s in screens}) == 6
    for path in [plan['union_plan'], union_plan['native_plan'], plan['union_completion'], plan['screens_plan']]: bind(bindings, path)
    # Native pins anchor the full inventory and written-input manifests used by
    # the independently verified union; no endpoint or denominator substitution.
    for path, digest in native_plan['pins'].items(): bind(bindings, path, digest)
    verify(bindings)
    return root / 'reference_measurement_dispositions.jsonl', pairs, models, inputs, receipt, bindings


def load_inventory_inputs(native_plan, plan, bindings):
    """Full static input preflight, also used before the union has finished."""
    inventory = Path(native_plan['inventory'])
    ir = json.loads((inventory / 'receipt.json').read_text())
    ip = json.loads(Path(native_plan['inventory_readback']).read_text())
    assert ir['status'] == 'complete_provisional_reference_comparison_inventory'
    assert ip['status'] == 'passed_full_reference_comparison_ledger_and_native_model_readback'
    assert ip['producer_receipt_sha256'] == sha(inventory / 'receipt.json')
    assert ip['distinct_model_pairs'] == ir['unique_distinct_model_pairs'] == plan['expected']['pairs']
    assert ip['unique_models'] == ir['all_reference_comparison_models'] == plan['expected']['models']
    artifacts(bindings, inventory, ir)
    bind(bindings, native_plan['inventory_readback'])
    models = {}
    with (inventory / 'models.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line); key = row['model_id'], row['version']
            assert key not in models and isinstance(row['version'], int) and row['length'] > 0
            models[key] = row
    assert len(models) == plan['expected']['models']
    pairs = {}
    with (inventory / 'model_pairs.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            ends = pair_ends(row)
            assert ends[0] != ends[1] and all(e in models for e in ends)
            assert row['pair_key'] == hashlib.sha256(json.dumps(sorted(ends), separators=(',', ':')).encode()).hexdigest()
            assert row['pair_key'] not in pairs
            pairs[row['pair_key']] = row
    assert len(pairs) == plan['expected']['pairs']
    # Scan both entire manifests, retaining only small per-input summaries.
    # Full protein lengths/sequence hashes are checked independently of masked lengths.
    inputs = {}; all_seen = set()
    for label, spec in native_plan['input_sources'].items():
        folder = Path(spec['inputs']); rp = folder / 'receipt.json'
        r = json.loads(rp.read_text()); p = json.loads(Path(spec['readback']).read_text())
        status = 'complete_duplication_alignment_input_materialization' if label == 'primary' else 'complete_additional_reference_alignment_input_materialization'
        proof_status = 'passed_full_duplication_alignment_input_readback' if label == 'primary' else 'passed_full_reference_alignment_input_readback'
        field = 'source_receipt_sha256' if label == 'primary' else 'producer_receipt_sha256'
        assert r['status'] == status and p['status'] == proof_status and p[field] == sha(rp)
        assert r['plan_sha256'] == sha(spec['input_plan'])
        assert json.loads(Path(spec['input_plan']).read_text())['output'] == str(folder)
        artifacts(bindings, folder, r)
        for path in [spec['input_plan'], spec['readback']]: bind(bindings, path)
        counts = Counter(); seen = set()
        with (folder / 'inputs.jsonl').open() as handle:
            for line in handle:
                row = json.loads(line); key = row['model_id'], row['version'], row['mask']
                assert key not in all_seen and key[2] in ['full', 'plddt70']
                assert row['status'] in ['ready', 'too_few_retained_residues', 'source_rejected']
                all_seen.add(key); seen.add(key); counts[key[2] + ':' + row['status']] += 1
                if key[:2] not in models: continue
                m = models[key[:2]]
                assert row['original_length'] == m['length'] and row['source_sha256'] == m['sha256']
                positions = row['original_positions']
                assert len(positions) == row['retained_residues'] == len(row['sequence'])
                assert positions == sorted(set(positions))
                assert all(isinstance(i, int) and 1 <= i <= m['length'] for i in positions)
                if key[2] == 'full' and row['status'] != 'source_rejected':
                    assert positions == list(range(1, m['length'] + 1))
                    assert hashlib.sha256(row['sequence'].encode()).hexdigest() == m['sequence_sha256']
                if row['status'] == 'ready': assert len(positions) >= 3
                elif row['status'] == 'too_few_retained_residues': assert len(positions) < 3
                inputs[key] = {k: row[k] for k in ['status', 'original_length', 'retained_residues']}
        assert len(seen) == r['input_dispositions'] == 2 * r['models'] and dict(counts) == r['counts']
        assert seen == {(m, v, mask) for m, v, _ in seen for mask in ['full', 'plddt70']}
    assert len(all_seen) == plan['expected']['all_input_dispositions']
    assert set(inputs) == {(*end, mask) for end in models for mask in ['full', 'plddt70']}
    return pairs, models, inputs


def groups(path):
    """Yield complete adjacent pair/mask groups from the sorted full union."""
    with Path(path).open() as handle:
        parsed = ((json.loads(line), hashlib.sha256(line.encode()).hexdigest()) for line in handle)
        for key, group in groupby(parsed, key=lambda r: (r[0]['pair_key'], r[0]['mask'])):
            yield key, list(group)
