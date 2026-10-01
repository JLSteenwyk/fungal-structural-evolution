"""Full expanded-background inventory, catalog partition and current input proof I/O."""
import csv
import json
from pathlib import Path
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

MASKS = ['full', 'plddt70']
PREFERENCE = ['primary', 'reference', 'background', 'legacy_reference']
SUMMARY_FIELDS = ['active_models', 'active_model_mask_states', 'full_pairs', 'new_pairs', 'pending_catalog_reuse_pairs',
                  'collection_model_counts', 'collection_input_states', 'overlapping_models', 'overlapping_model_mask_states',
                  'active_input_counts', 'directed_native_input_counts', 'actual_ready_PDB_files', 'actual_raw_coordinate_models',
                  'raw_coordinate_bytes', 'ready_PDB_bytes']


def load_sources(plan, plan_path):
    bindings = dict(plan['pins']); bind(bindings, plan_path)
    completion = json.loads(Path(plan['background_input_completion']).read_text())
    assert completion['status'] == 'complete_verified_full_expanded_background_native_coordinates_and_PDB_inputs' and len(completion['services']) == 4
    for path, digest in completion['source_hashes'].items(): bind(bindings, path, digest)
    bind(bindings, plan['background_input_completion'])
    legacy = json.loads(Path(plan['legacy_input_completion']).read_text())
    assert legacy['status'] == 'complete_verified_full_legacy_reference_written_inputs' and len(legacy['services']) == 2
    for path, digest in legacy['source_hashes'].items(): bind(bindings, path, digest)
    bind(bindings, plan['legacy_input_completion'])
    inventory = Path(plan['inventory']); rp = inventory / 'receipt.json'; r = json.loads(rp.read_text())
    proof = json.loads(Path(plan['inventory_readback']).read_text())
    assert r['status'] == 'complete_background_measurement_inventory_pending_readback'
    assert proof['status'] == 'passed_full_background_measurement_inventory_readback' and proof['producer_receipt_sha256'] == sha(rp)
    assert r['active_models'] == proof['active_models'] == plan['expected']['active_models']
    assert r['distinct_eligible_model_pairs'] == proof['distinct_eligible_model_pairs'] == plan['expected']['full_pairs']
    bind(bindings, rp); bind(bindings, plan['inventory_readback'])
    for name, digest in r['artifacts'].items(): bind(bindings, inventory / name, digest)
    active = {}
    with (inventory / 'active_models.jsonl').open() as handle:
        for line in handle:
            m = json.loads(line); key = m['model_id'], m['version']; assert key not in active; active[key] = m
    assert len(active) == plan['expected']['active_models']
    with (inventory / 'model_pairs.tsv').open() as handle: pairs = list(csv.DictReader(handle, delimiter='\t'))
    assert len(pairs) == plan['expected']['full_pairs']
    root = Path(plan['reuse_candidates']); rr = json.loads((root / 'receipt.json').read_text()); ra = json.loads(Path(plan['reuse_readback']).read_text())
    assert rr['status'] == 'complete_full_pair_union_catalog_reuse_screen_not_authorization'
    assert ra['status'] == 'passed_full_pair_union_reuse_candidate_sql_readback' and ra['producer_receipt_sha256'] == sha(root / 'receipt.json')
    bind(bindings, root / 'receipt.json'); bind(bindings, plan['reuse_readback'])
    for name, digest in rr['artifacts'].items(): bind(bindings, root / name, digest)
    for name in ['receipt.json', 'models.jsonl', 'model_pairs.tsv']: assert rr['source_hashes'][str(inventory / name)] == sha(inventory / name)
    candidates = {}
    with (root / 'pair_reuse_candidates.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if 'background_expanded' in json.loads(row['new_sources']): assert row['pair_key'] not in candidates; candidates[row['pair_key']] = row
    assert len(candidates) == len(pairs)
    specifications = {}
    assert list(plan['input_sources']) == PREFERENCE
    for label, spec in plan['input_sources'].items():
        folder = Path(spec['inputs']); irp = folder / 'receipt.json'; ir = json.loads(irp.read_text()); config = json.loads(Path(spec['input_plan']).read_text())
        assert ir['status'] == spec['status'] and ir['plan_sha256'] == sha(spec['input_plan']) and Path(config['output']).resolve() == folder.resolve()
        source = Path(spec['model_inventory']); mrp = source / 'receipt.json'; mr = json.loads(mrp.read_text())
        assert ir[spec['inventory_receipt_field']] == sha(mrp)
        bind(bindings, irp); bind(bindings, mrp); bind(bindings, source / spec['model_file'], mr['artifacts'][spec['model_file']])
        bind(bindings, folder / 'inputs.jsonl', ir['artifacts']['inputs.jsonl']); bind(bindings, spec['input_plan'])
        audit = json.loads(Path(spec['readback']).read_text())
        assert audit['status'] == spec['readback_status'] and audit[spec['readback_receipt_field']] == sha(irp)
        assert audit['models'] == ir['models'] and audit['input_dispositions'] == ir['input_dispositions'] and audit['counts'] == ir['counts']
        bind(bindings, spec['readback'])
        coordinate = Path(config['coordinates']) / 'receipt.json'; native = Path(config['readback']) / 'receipt.json'; nr = json.loads(native.read_text())
        assert ir['coordinate_receipt_sha256'] == nr['producer_receipt_sha256'] == sha(coordinate)
        assert ir['readback_receipt_sha256'] == sha(native) and nr['status'] == spec['coordinate_readback_status']
        bind(bindings, coordinate); bind(bindings, native)
        specifications[label] = dict(spec=spec, receipt=ir, model_path=source / spec['model_file'], manifest=folder / 'inputs.jsonl')
    verify(bindings)
    return active, pairs, candidates, specifications, bindings
