"""Closed full background measurements, full model lengths and fixed target screens."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from background_measurement_union_sources import closed_source, SUMMARY_FIELDS as UNION_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

ORDER_FIELDS = ['selected_source', 'source_order', 'source_checkpoint_sha256', 'native_status', 'numerical_usable', 'numerical_exclusion_reasons', 'aligned_length', 'original_coverage_a', 'original_coverage_b']
SUMMARY_FIELDS = ['source_inventory_models', 'pair_models', 'pairs', 'pair_mask_rows', 'directed_source_states', 'screen_decisions', 'screens', 'pass_counts', 'both_masks_pass_counts', 'exclusion_counts', 'target_pairs', 'target_pass_counts', 'target_both_masks_pass_counts']


def load(plan, plan_path):
    bindings = dict(plan['pins']); bind(bindings, plan_path)
    c = closed_source(plan['measurement_completion'], 'complete_verified_full_background_measurement_union', 'complete_verified_full_background_measurement_union_archive', 2, bindings)
    assert c['scientific_eligibility'] is False and c['full_pairs'] == plan['expected']['pairs'] and c['directed_dispositions'] == 4 * c['full_pairs']
    assert bindings[c['producer_receipt']] == c['producer_receipt_sha256'] and bindings[c['independent_readback']] == c['independent_readback_sha256']
    assert c['source_plan'] == plan['measurement_plan'] and c['source_plan_sha256'] == sha(plan['measurement_plan'])
    up = json.loads(Path(plan['measurement_plan']).read_text()); root = Path(up['output'])
    rp, ap = Path(c['producer_receipt']), Path(c['independent_readback']); r, a = [json.loads(p.read_text()) for p in [rp, ap]]
    assert rp.parent == root and r['status'] == 'complete_full_background_measurement_union_pending_independent_readback' and a['status'] == 'passed_full_background_measurement_union_sql_readback'
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['measurement_plan']) and a['producer_receipt_sha256'] == sha(rp)
    assert all(r[key] == a[key] == c[key] for key in UNION_FIELDS)
    for name, digest in r['artifacts'].items(): bind(bindings, root / name, digest)
    model_path = Path(plan['models']); assert str(model_path) in bindings
    models = {}
    with gzip.open(model_path, 'rt') as handle:
        for line in handle:
            row = json.loads(line); key = row['model_id'], row['version']; assert key not in models and isinstance(row['length'], int) and row['length'] > 0
            assert math.isfinite(row['mean_ca_plddt']) and 0 <= row['mean_ca_plddt'] <= 100 and math.isfinite(row['fraction_ca_plddt_below50']) and 0 <= row['fraction_ca_plddt_below50'] <= 1
            models[key] = row
    assert len(models) == plan['expected']['models']
    ledger = root / 'full_background_work_partition.tsv'; pairs = {}
    with ledger.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            ends = [(row['model_a'], int(row['version_a'])), (row['model_b'], int(row['version_b']))]
            assert row['pair_key'] not in pairs and ends == sorted(ends) and ends[0] != ends[1] and all(end in models for end in ends)
            assert row['pair_key'] == hashlib.sha256(json.dumps(ends, separators=(',', ':')).encode()).hexdigest()
            pairs[row['pair_key']] = ends
    used_models = {end for ends in pairs.values() for end in ends}
    assert len(pairs) == c['full_pairs'] and used_models <= set(models) and len(used_models) == plan['expected']['pair_models']
    tc = Path(plan['target_completion']); target = json.loads(tc.read_text()); bind(bindings, tc)
    assert target['status'] == 'complete_verified_full_expanded_pair_event_and_taxon_coverage_screens' and len(target['services']) == 2 and target['scientific_eligibility'] is False
    for path, digest in target['source_hashes'].items(): bind(bindings, path, digest)
    tp = json.loads(Path(plan['target_plan']).read_text()); tr = Path(tp['output']) / 'receipt.json'; t = json.loads(tr.read_text())
    bind(bindings, plan['target_plan']); bind(bindings, plan['screens_plan']); bind(bindings, plan['measurement_plan'])
    assert str(tr) in bindings and t['status'] == 'complete_full_expanded_duplication_coverage_pending_independent_readback' and t['plan_sha256'] == sha(plan['target_plan'])
    assert t['screens'] == tp['screens'] == plan['screens'] == json.loads(Path(plan['screens_plan']).read_text())['screens']
    assert t['pairs'] == target['summary']['pairs'] == plan['expected']['target_pairs']
    for name in ['pair_mask_rows', 'pair_pass_counts', 'pair_both_masks_pass_counts']: assert t[name] == target['summary'][name]
    assert t['pair_mask_rows'] == 2 * t['pairs']
    verify(bindings)
    return dict(models=models, pair_models=len(used_models), pairs=pairs, states=root / 'background_measurement_dispositions.jsonl.gz', target=t), bindings


def fields(screens):
    return (['pair_key', 'mask', 'model_a', 'version_a', 'model_b', 'version_b', 'length_a', 'length_b', 'mean_ca_plddt_a', 'mean_ca_plddt_b', 'fraction_ca_plddt_below50_a', 'fraction_ca_plddt_below50_b']
            + [f'order{order}_{name}' for order in [0, 1] for name in ORDER_FIELDS]
            + [s['id'] + suffix for s in screens for suffix in ['_pass', '_exclusions']])
