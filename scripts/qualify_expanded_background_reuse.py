#!/usr/bin/env python3
"""Qualify every old background candidate with exact physical inputs and native results."""
import argparse
import gzip
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from duplication_alignment_numeric_diagnostic import check_alignment
from duplication_alignment_numeric_readback import load_pdb
from expanded_background_reuse_sources import load_sources
from qualify_reference_alignment_reuse import endpoints, projection, digest, inspect_checkpoint, numeric_exclusions
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text())
    full, selected, current, models, old, bindings, source_root = load_sources(plan, plan_path)
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    ledger = 'full_background_work_partition.tsv'; (out / ledger).write_bytes((source_root / ledger).read_bytes())
    checks = {}; checked_files = set(); totals, numeric_counts, owners = Counter(), Counter(), Counter(); numeric_checked = 0
    @lru_cache(maxsize=256)
    def coordinates(key): return load_pdb(current[key])
    def actual(path, expected):
        if (path, expected) not in checked_files: assert sha(path) == expected; checked_files.add((path, expected))
        bind(bindings, path, expected)
    @lru_cache(maxsize=None)
    def model_check(owner, key, mask):
        a, b = current[(*key, mask)], old[owner]['inputs'][(*key, mask)]; model, previous = models[key], old[owner]['models'][key]
        compatible = all(model[name] == previous[name] for name in ['sha256', 'sequence_sha256', 'length']) and projection(a) == projection(b)
        for row in [model, previous]: actual(row['path'], row['sha256'])
        for row, raw in [(a, model), (b, previous)]:
            assert row['source_sha256'] == raw['sha256']
            if row['status'] == 'ready': actual(row['path'], row['sha256'])
        checks[(owner, *key, mask)] = dict(source=owner, model_id=key[0], version=key[1], mask=mask, compatible=compatible,
            current_projection_sha256=digest(projection(a)), old_projection_sha256=digest(projection(b)), current_raw_path=model['path'], old_raw_path=previous['path'],
            current_raw_sha256=model['sha256'], old_raw_sha256=previous['sha256'], current_pdb_path=a.get('path'), old_pdb_path=b.get('path'))
        return compatible
    with gzip.open(out / 'background_reuse_dispositions.jsonl.gz', 'wt') as handle:
        for index, (pair, original) in enumerate(sorted(full.items()), 1):
            ends = endpoints(original); owner = selected[pair]
            for mask in ['full', 'plddt70']:
                compatible = bool(owner) and all([model_check(owner, end, mask) for end in ends])
                for order in [0, 1]:
                    row = dict(pair_key=pair, model_a=ends[0][0], version_a=ends[0][1], model_b=ends[1][0], version_b=ends[1][1], mask=mask, order=order,
                        matching_old_sources=json.loads(original['matching_old_sources']), selected_source=owner, numerical_usable=False, source_checkpoint=None, source_checkpoint_sha256=None,
                        source_order=None, source_native_status=None, source_numeric=None, source_geometry=None, numerical_exclusion_reasons=[])
                    if not owner: row['reuse_status'] = 'new_native_measurement_pending'
                    elif not compatible: row['reuse_status'] = 'incompatible_inputs_require_new_measurement'
                    else:
                        source = old[owner]; desired = ends if order == 0 else ends[::-1]
                        mappings = [source['pairs'][pair], source['pairs'][pair][::-1]]; assert mappings.count(desired) == 1; source_order = mappings.index(desired)
                        path, h, native, numeric, geometry = inspect_checkpoint(source, pair, mask, source_order, bindings)
                        if native['status'] == 'aligned':
                            reconstructed = check_alignment(native, *[coordinates((*end, mask)) for end in desired])
                            assert set(numeric) == {'pair_key', 'mask', 'order', *reconstructed}
                            for name, value in reconstructed.items():
                                assert numeric[name] == value if isinstance(value, str) else float(numeric[name]) == float(value), name
                            numeric_checked += 1
                        reasons = numeric_exclusions(native, numeric, geometry)
                        row.update(reuse_status='verified_identical_input_checkpoint_and_retained_disposition', source_checkpoint=str(path), source_checkpoint_sha256=h, source_order=source_order,
                                   source_native_status=native['status'], source_numeric=numeric, source_geometry=geometry, numerical_exclusion_reasons=reasons, numerical_usable=not reasons)
                        numeric_counts['usable' if not reasons else ';'.join(reasons)] += 1
                    totals[row['reuse_status']] += 1; owners[owner or 'no_old_source'] += 1; handle.write(json.dumps(row, separators=(',', ':')) + '\n')
            if index % 1000 == 0:
                (out / 'state.json').write_text(json.dumps(dict(stage='qualifying_full_background_reuse', pairs=index, total_pairs=len(full), counts=dict(totals))) + '\n')
                print('Full background actual reuse pairs', index, '/', len(full), flush=True)
    with gzip.open(out / 'model_input_identity_checks.jsonl.gz', 'wt') as handle:
        for _, row in sorted(checks.items()): handle.write(json.dumps(row, separators=(',', ':')) + '\n')
    verify(bindings)
    result = dict(status='complete_full_expanded_background_reuse_qualification_pending_independent_readback', plan_sha256=sha(plan_path), full_background_pairs=len(full), directed_dispositions=4 * len(full),
                  counts=dict(totals), numerical_counts=dict(numeric_counts), selected_source_dispositions=dict(owners), unique_model_mask_source_checks=len(checks), numerically_reconstructed_reuse_alignments=numeric_checked,
                  source_hashes=bindings, artifacts={name: sha(out / name) for name in ['background_reuse_dispositions.jsonl.gz', 'model_input_identity_checks.jsonl.gz', ledger]}, scientific_eligibility=False, scope=plan['scope'])
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); run(parser.parse_args().plan)
