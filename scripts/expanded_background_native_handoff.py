"""Consume the closed full input/work union; preserve pending actual result reuse."""
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def load_handoff(plan, include_positions=False):
    bindings = dict(plan['pins']); cp = plan['handoff_completion']; c = json.loads(Path(cp).read_text()); bind(bindings, cp)
    assert c['status'] == 'complete_verified_full_expanded_background_native_input_handoff' and c['exact_process_journals_checked'] == 2
    bind(bindings, c['full_hash_archive'], c['full_hash_archive_sha256']); archive = json.loads(Path(c['full_hash_archive']).read_text())
    assert archive['status'] == 'complete_verified_full_expanded_background_native_input_handoff_archive' and len(archive['services']) == 2
    assert len(archive['source_hashes']) == c['bound_source_hashes']
    for path, digest in archive['source_hashes'].items(): bind(bindings, path, digest)
    rp, ap = Path(c['producer_receipt']), Path(c['independent_readback'])
    bind(bindings, rp, c['producer_receipt_sha256']); bind(bindings, ap, c['independent_readback_sha256'])
    r, a = [json.loads(p.read_text()) for p in [rp, ap]]
    assert r['status'] == 'complete_full_expanded_background_native_input_handoff_pending_independent_readback'
    assert a['status'] == 'passed_full_expanded_background_native_input_handoff_sql_readback' and a['producer_receipt_sha256'] == sha(rp)
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['handoff_plan'])
    for key in ['active_models', 'active_model_mask_states', 'full_pairs', 'new_pairs', 'pending_catalog_reuse_pairs']:
        assert r[key] == a[key] == c[key] == plan['expected'][key]
    bind(bindings, plan['handoff_plan']); root = rp.parent
    for name, digest in r['artifacts'].items(): bind(bindings, root / name, digest)
    inputs = {}; origins = {}; statuses = Counter()
    with gzip.open(root / 'active_inputs.jsonl.gz', 'rt') as handle:
        for line in handle:
            row = json.loads(line); key = row['model_id'], row['version'], row['mask']; assert key not in inputs and key[2] in ['full', 'plddt70']
            inputs[key] = row; statuses[key[2] + ':' + row['status']] += 1
            if include_positions:
                origin = row['collection_sources'][0]; source_key = origin['source_manifest'], origin['source_row_number']; assert source_key not in origins; origins[source_key] = key
    assert len(inputs) == c['active_model_mask_states'] == 2 * c['active_models'] and dict(statuses) == c['active_input_counts']
    if include_positions:
        config = json.loads(Path(plan['handoff_plan']).read_text()); found = set()
        for spec in config['input_sources'].values():
            manifest = str(Path(spec['inputs']) / 'inputs.jsonl'); assert manifest in bindings
            with Path(manifest).open() as handle:
                for ordinal, line in enumerate(handle, 1):
                    source_key = manifest, ordinal
                    if source_key not in origins: continue
                    key = origins[source_key]; row = json.loads(line); selected = inputs[key]
                    assert (row['model_id'], row['version'], row['mask']) == key
                    fingerprint = hashlib.sha256(json.dumps(row, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
                    assert fingerprint == selected['collection_sources'][0]['source_row_sha256']
                    semantic = {k: v for k, v in row.items() if k not in ['path', 'coordinate_shard']}
                    assert hashlib.sha256(json.dumps(semantic, sort_keys=True, separators=(',', ':')).encode()).hexdigest() == selected['semantic_sha256']
                    for name in ['status', 'source_sha256', 'path', 'sha256', 'sequence', 'retained_residues', 'original_length']:
                        assert row.get(name) == selected.get(name)
                    inputs[key]['original_positions'] = row.get('original_positions', []); found.add(source_key)
        assert found == set(origins)
    partition = []; new = []; old = 0; keys = set()
    with (root / 'full_background_work_partition.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['pair_key']; assert key not in keys; keys.add(key)
            ends = [(row['model_a'], int(row['version_a'])), (row['model_b'], int(row['version_b']))]
            assert ends == sorted(ends) and ends[0] != ends[1] and hashlib.sha256(json.dumps(ends, separators=(',', ':')).encode()).hexdigest() == key
            assert all((*end, mask) in inputs for end in ends for mask in ['full', 'plddt70'])
            if row['measurement_disposition'] == 'native_measurement_required': assert not json.loads(row['matching_old_sources']); new.append(row)
            else: assert row['measurement_disposition'] == 'pending_actual_input_result_and_numeric_reuse_checks' and json.loads(row['matching_old_sources']); old += 1
            partition.append(row)
    assert len(partition) == c['full_pairs'] and len(new) == c['new_pairs'] and old == c['pending_catalog_reuse_pairs']
    verify(bindings)
    return inputs, new, bindings, r['artifacts']['active_inputs.jsonl.gz'], partition, root
