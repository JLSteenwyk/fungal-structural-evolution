"""Full current input union and original reference/background native source contracts."""
import csv
import gzip
import json
from collections import Counter
from pathlib import Path
from background_alignment_handoff import load_handoff as old_background_handoff
from expanded_background_native_handoff import load_handoff
from qualify_reference_alignment_reuse import load_source, table_subset, endpoints
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def old_background(spec, wanted, bindings):
    sp = Path(spec['alignment_plan']); plan = json.loads(sp.read_text()); root = Path(plan['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['status'] == spec['alignment_status'] and receipt['plan_sha256'] == sha(sp)
    compact, measured, source_bindings, bundle = old_background_handoff(plan)
    assert source_bindings == receipt['input_bindings'] and bundle == receipt['input_bundle_sha256']
    for path, digest in source_bindings.items(): bind(bindings, path, digest)
    pairs = {row['pair_key']: endpoints(row) for row in measured}; assert len(pairs) == receipt['distinct_model_pairs']
    assert wanted <= set(pairs) and receipt['directed_dispositions'] == 4 * len(pairs)
    wanted_models = {end for pair in wanted for end in pairs[pair]}; models = {}
    mp = Path(plan['inventory']) / 'active_models.jsonl'; bind(bindings, mp)
    with mp.open() as handle:
        for line in handle:
            row = json.loads(line); key = row['model_id'], row['version']
            if key in wanted_models: assert key not in models; models[key] = row
    assert set(models) == wanted_models
    inputs = {}
    for config in plan['input_sources'].values():
        manifest = Path(config['inputs']) / 'inputs.jsonl'; assert str(manifest) in source_bindings
        with manifest.open() as handle:
            for line in handle:
                row = json.loads(line); key = row['model_id'], row['version'], row['mask']
                if key[:2] in wanted_models:
                    assert key not in inputs and row['source_sha256'] == models[key[:2]]['sha256']
                    assert {name: row[name] for name in ['status', 'path', 'sha256', 'sequence'] if name in row} == compact[key]
                    inputs[key] = row
    del compact
    assert set(inputs) == {(*model, mask) for model in wanted_models for mask in ['full', 'plddt70']}
    proofs = {}; full_keys = set(); counts = Counter(); manifest = root / 'checkpoint_manifest.tsv'
    bind(bindings, root / 'receipt.json')
    for name, digest in receipt['artifacts'].items(): bind(bindings, root / name, digest)
    with manifest.open() as handle:
        for item in csv.DictReader(handle, delimiter='\t'):
            pair, mask, order = Path(item['path']).stem.rsplit('-', 2); key = pair, mask, int(order)
            assert key not in full_keys and pair in pairs and mask in ['full', 'plddt70'] and key[2] in [0, 1]
            assert item['path'] == f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json'
            full_keys.add(key); counts[mask + ':' + item['status']] += 1
            if pair in wanted: proofs[key] = item
    assert full_keys == {(p, m, o) for p in pairs for m in ['full', 'plddt70'] for o in [0, 1]} and dict(counts) == receipt['counts']
    dp, gp = Path(spec['diagnostic']), Path(spec['geometry']); dr, gr = [json.loads((p / 'receipt.json').read_text()) for p in [dp, gp]]
    numeric = table_subset(dp / 'numeric_readback.tsv', wanted, None, dr['numerically_checked_alignments'])
    geometry = table_subset(gp / 'alignment_geometry.tsv', wanted, None, gr['alignments']); assert set(numeric) == set(geometry)
    bind(bindings, sp)
    return dict(plan=plan, plan_sha256=sha(sp), root=root, pairs=pairs, models=models, inputs=inputs, proofs=proofs, numeric=numeric, geometry=geometry, bundle=bundle)


def proof_lineage(spec, source, bindings):
    sp, dp, gp, qp = [Path(spec[key]) for key in ['alignment_plan', 'diagnostic_plan', 'geometry_plan', 'geometry_readback_plan']]
    configs = [json.loads(p.read_text()) for p in [sp, dp, gp, qp]]
    assert configs[1]['source_plan'] == str(sp) and configs[2]['diagnostic_plan'] == str(dp) and configs[3]['source_plan'] == str(gp)
    roots = [Path(config['output']) for config in configs[:3]]; assert roots[0] == source['root']
    assert str(roots[1]) == spec['diagnostic'] and str(roots[2]) == spec['geometry'] and configs[3]['output'] == spec['geometry_readback']
    native, diagnostic, geometry = [json.loads((root / 'receipt.json').read_text()) for root in roots]; reader = json.loads(Path(spec['geometry_readback']).read_text())
    for root, record, config in zip(roots, [native, diagnostic, geometry], [sp, dp, gp]):
        assert record['plan_sha256'] == sha(config); bind(bindings, root / 'receipt.json')
        for name, digest in record['artifacts'].items(): bind(bindings, root / name, digest)
    assert diagnostic['status'] == spec['diagnostic_status'] and geometry['status'] == spec['geometry_status'] and reader['status'] == spec['geometry_readback_status']
    assert diagnostic['producer_receipt_sha256'] == sha(roots[0] / 'receipt.json') and geometry['diagnostic_receipt_sha256'] == sha(roots[1] / 'receipt.json')
    assert reader['producer_receipt_sha256'] == sha(roots[2] / 'receipt.json') and reader['plan_sha256'] == sha(qp)
    assert diagnostic['directed_dispositions'] == native['directed_dispositions'] and diagnostic['counts'] == native['counts']
    assert diagnostic['numerically_checked_alignments'] == geometry['alignments'] == reader['alignments_checked']
    assert sum(value for key, value in native['counts'].items() if key.endswith(':aligned')) == geometry['alignments'] and reader['counts'] == geometry['counts']
    aligned = {key for key, item in source['proofs'].items() if item['status'] == 'aligned'}
    assert set(source['numeric']) == set(source['geometry']) == aligned
    for p in [sp, dp, gp, qp, Path(spec['geometry_readback']), *[Path(p) for p in spec['launches']]]: bind(bindings, p)


def load_sources(plan, plan_path):
    native_path = Path(plan['target_alignment_plan']); target = json.loads(native_path.read_text())
    current, new, bindings, _, partition, root = load_handoff(target, include_positions=True)
    bind(bindings, native_path); bind(bindings, plan_path)
    for path, digest in plan['pins'].items(): bind(bindings, path, digest)
    models = {}
    with gzip.open(root / 'active_models.jsonl.gz', 'rt') as handle:
        for line in handle:
            row = json.loads(line); key = row['model_id'], row['version']; assert key not in models; models[key] = row
    assert {key[:2] for key in current} == set(models) and len(models) == plan['expected']['active_models']
    full = {row['pair_key']: row for row in partition}; selected = {}; source_counts = Counter()
    for pair, row in full.items():
        matches = json.loads(row['matching_old_sources']); assert not (set(matches) - {'reference_old', 'background_old'})
        choice = next((name for name in ['reference_old', 'background_old'] if name in matches), '')
        assert bool(choice) == (row['measurement_disposition'] == 'pending_actual_input_result_and_numeric_reuse_checks')
        selected[pair] = choice; source_counts[choice or 'no_old_source'] += 1
    assert len(full) == plan['expected']['full_pairs'] and source_counts['no_old_source'] == len(new) == plan['expected']['new_pairs']
    assert {k: v for k, v in source_counts.items() if k != 'no_old_source'} == plan['expected']['selected_source_pairs']
    sources = {}
    for label, spec in plan['sources'].items():
        wanted = {pair for pair, owner in selected.items() if owner == label}
        sources[label] = old_background(spec, wanted, bindings) if spec['kind'] == 'background' else load_source(spec, wanted, bindings)
        proof_lineage(spec, sources[label], bindings)
        old = sources[label]['plan']; h = sha(target['usalign'])
        assert old['options'] == target['options'] and old['per_pair_timeout_seconds'] == target['per_pair_timeout_seconds'] and sha(old['usalign']) == old['pins'][old['usalign']] == h
        for pair in wanted: assert sorted(sources[label]['pairs'][pair]) == endpoints(full[pair])
    verify(bindings)
    return full, selected, current, models, sources, bindings, root
