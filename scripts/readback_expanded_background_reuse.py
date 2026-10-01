#!/usr/bin/env python3
"""Rebuild the full reuse grid and every accepted original mapping with quaternion fits."""
import argparse
import gzip
import json
import math
from collections import Counter
from functools import lru_cache
from pathlib import Path
from duplication_alignment_numeric_readback import load_pdb
from expanded_background_reuse_sources import load_sources
from readback_reference_alignment_reuse import signature, identity_hash, native_fields
from readback_expanded_background_measurements import numeric_fields
from readback_background_alignment_geometry import verify_row
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path, output):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text()); root = Path(plan['output']); rp = root / 'receipt.json'; r = json.loads(rp.read_text()); rh = sha(rp)
    assert r['status'] == 'complete_full_expanded_background_reuse_qualification_pending_independent_readback' and r['plan_sha256'] == sha(plan_path) and r['scientific_eligibility'] is False
    bindings = {}; bind(bindings, plan_path)
    for path, digest in r['source_hashes'].items(): bind(bindings, path, digest)
    bind(bindings, rp, rh)
    for name, digest in r['artifacts'].items(): bind(bindings, root / name, digest)
    verify(bindings)
    full, loader_choices, current, models, sources, loaded, source_root = load_sources(plan, plan_path)
    assert all(bindings.get(path) == digest for path, digest in loaded.items())
    assert (root / 'full_background_work_partition.tsv').read_bytes() == (source_root / 'full_background_work_partition.tsv').read_bytes()
    chosen = {}
    for pair, row in full.items():
        candidates = json.loads(row['matching_old_sources'])
        chosen[pair] = 'reference_old' if 'reference_old' in candidates else 'background_old' if 'background_old' in candidates else ''
    assert chosen == loader_choices
    model_checks = {}; reconstructed_bindings = dict(loaded)
    def actual(path, expected):
        assert bindings.get(path) == expected
        bind(reconstructed_bindings, path, expected)
    for pair, original in full.items():
        owner = chosen[pair]
        if not owner: continue
        for end in [(original['model_a'], int(original['version_a'])), (original['model_b'], int(original['version_b']))]:
            for mask in ['full', 'plddt70']:
                key = owner, *end, mask
                if key in model_checks: continue
                a, b = current[(*end, mask)], sources[owner]['inputs'][(*end, mask)]; model, old = models[end], sources[owner]['models'][end]
                equal = (model['sha256'], model['sequence_sha256'], model['length']) == (old['sha256'], old['sequence_sha256'], old['length']) and signature(a) == signature(b)
                for value in [model, old]: actual(value['path'], value['sha256'])
                for value, raw in [(a, model), (b, old)]:
                    assert value['source_sha256'] == raw['sha256']
                    if value['status'] == 'ready': actual(value['path'], value['sha256'])
                model_checks[key] = dict(source=owner, model_id=end[0], version=end[1], mask=mask, compatible=equal,
                    current_projection_sha256=identity_hash(a), old_projection_sha256=identity_hash(b), current_raw_path=model['path'], old_raw_path=old['path'],
                    current_raw_sha256=model['sha256'], old_raw_sha256=old['sha256'], current_pdb_path=a.get('path'), old_pdb_path=b.get('path'))
    with gzip.open(root / 'model_input_identity_checks.jsonl.gz', 'rt') as handle:
        for _, expected in sorted(model_checks.items()): assert json.loads(next(handle)) == expected
        assert next(handle, None) is None
    counts, numerical, owners = Counter(), Counter(), Counter(); checked = numeric_checked = near = 0; maximum = maximum_curvature = 0.
    @lru_cache(maxsize=256)
    def coordinates(key): return load_pdb(current[key])
    with gzip.open(root / 'background_reuse_dispositions.jsonl.gz', 'rt') as handle:
        for pair, original in sorted(full.items()):
            ends = [(original['model_a'], int(original['version_a'])), (original['model_b'], int(original['version_b']))]; owner = chosen[pair]
            for mask in ['full', 'plddt70']:
                compatible = bool(owner) and all(model_checks[(owner, *end, mask)]['compatible'] for end in ends)
                for order in [0, 1]:
                    row = json.loads(next(handle))
                    expected = dict(pair_key=pair, model_a=ends[0][0], version_a=ends[0][1], model_b=ends[1][0], version_b=ends[1][1], mask=mask, order=order,
                        matching_old_sources=json.loads(original['matching_old_sources']), selected_source=owner, numerical_usable=False, source_checkpoint=None, source_checkpoint_sha256=None,
                        source_order=None, source_native_status=None, source_numeric=None, source_geometry=None, numerical_exclusion_reasons=[])
                    if not owner: expected['reuse_status'] = 'new_native_measurement_pending'
                    elif not compatible: expected['reuse_status'] = 'incompatible_inputs_require_new_measurement'
                    else:
                        source = sources[owner]; desired = ends if order == 0 else list(reversed(ends))
                        mappings = [source['pairs'][pair], list(reversed(source['pairs'][pair]))]; assert mappings.count(desired) == 1; old_order = mappings.index(desired)
                        expected.update(native_fields(source, pair, mask, old_order, desired)); expected['reuse_status'] = 'verified_identical_input_checkpoint_and_retained_disposition'
                        actual(expected['source_checkpoint'], expected['source_checkpoint_sha256'])
                        if expected['source_native_status'] == 'aligned':
                            record = json.loads(Path(expected['source_checkpoint']).read_text()); numeric, x, y = numeric_fields(record, [coordinates((*end, mask)) for end in desired])
                            original_numeric = expected['source_numeric']; assert set(original_numeric) == {'pair_key', 'mask', 'order', *numeric}
                            for name, value in numeric.items():
                                if isinstance(value, str): assert original_numeric[name] == value
                                else: assert math.isclose(float(original_numeric[name]), value, rel_tol=1e-9, abs_tol=1e-9), name
                            maximum = max(maximum, abs(float(original_numeric['rmsd_recomputed']) - numeric['rmsd_recomputed']))
                            err, boundary = verify_row(expected['source_geometry'], x, y); maximum_curvature = max(maximum_curvature, err); near += boundary; numeric_checked += 1
                            if mask == 'plddt70': assert numeric['joint_plddt70_pairs'] == numeric['aligned_length']
                        numerical['usable' if expected['numerical_usable'] else ';'.join(expected['numerical_exclusion_reasons'])] += 1
                    assert row == expected, pair
                    counts[row['reuse_status']] += 1; owners[owner or 'no_old_source'] += 1; checked += 1
                    if checked % 10000 == 0: print('Independent full background actual reuse states', checked, flush=True)
        assert next(handle, None) is None
    summary = dict(full_background_pairs=len(full), directed_dispositions=checked, counts=dict(counts), numerical_counts=dict(numerical), selected_source_dispositions=dict(owners),
                   unique_model_mask_source_checks=len(model_checks), numerically_reconstructed_reuse_alignments=numeric_checked)
    assert checked == 4 * len(full) and all(r[key] == value for key, value in summary.items())
    assert reconstructed_bindings == r['source_hashes']; verify(bindings)
    result = dict(status='passed_full_expanded_background_input_checkpoint_numeric_quaternion_reuse_readback', plan_sha256=sha(plan_path), producer_receipt_sha256=rh,
                  **summary, maximum_absolute_reused_rmsd_difference=maximum, maximum_scaled_reused_quaternion_curvature_error=maximum_curvature, near_zero_reused_quaternion_gaps=near,
                  source_hashes=bindings, scientific_eligibility=False, scope=plan['scope'])
    with Path(output).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); run(args.plan, args.output)
