#!/usr/bin/env python3
"""Reconstruct every new background disposition, numeric mapping and rigid geometry."""
import argparse
import gzip
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
import numpy as np
from assess_background_alignment_geometry import geometry
from duplication_alignment_numeric_diagnostic import check_alignment
from duplication_alignment_numeric_readback import load_pdb
from expanded_background_measurement_sources import load_native, inspect_native
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text())
    source, root, receipt, inputs, pairs, checkpoints, bindings, bundle = load_native(plan)
    native_plan_hash = sha(plan['native_plan'])
    bind(bindings, plan_path); out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    ledger = 'full_background_work_partition.tsv'; (out / ledger).write_bytes((root / ledger).read_bytes())
    native_counts, numeric_counts, geometry_counts, usable_counts = [Counter() for _ in range(4)]
    checked = 0; maximum = 0.
    @lru_cache(maxsize=256)
    def coordinates(key): return load_pdb(inputs[key])
    with gzip.open(out / 'disposition_measurements.jsonl.gz', 'wt') as handle:
        for index, (key, item) in enumerate(checkpoints.items(), 1):
            pair, mask, order = key; ends = pairs[pair] if order == 0 else pairs[pair][::-1]
            rows = [inputs[(*end, mask)] for end in ends]
            record = inspect_native(source, root, item, key, ends, rows, bundle)
            assert record['plan_sha256'] == native_plan_hash
            result = dict(pair_key=pair, mask=mask, order=order, model_left=ends[0][0], version_left=ends[0][1], model_right=ends[1][0], version_right=ends[1][1],
                          checkpoint_path=str(root / item['path']), checkpoint_sha256=item['sha256'], native_status=record['status'], native_elapsed_seconds=record['elapsed_seconds'],
                          numerical=None, geometry=None, numerical_usable=False, numerical_exclusion_reasons=[])
            if record['status'] == 'aligned':
                coords = [coordinates((*end, mask)) for end in ends]; numeric = check_alignment(record, *coords)
                strings = [record['metrics'][name] for name in ['alignment_left', 'alignment_right']]
                indexes = [[], []]; counters = [0, 0]
                for chars in zip(*strings):
                    if all(char != '-' for char in chars):
                        for side in range(2): indexes[side].append(counters[side])
                    for side in range(2): counters[side] += chars[side] != '-'
                g = geometry(*[c[1][ix] for c, ix in zip(coords, indexes)]); reasons = []
                if numeric['rmsd_status'] != 'within_printed_rounding': reasons.append('rmsd_discrepancy')
                if numeric['aligned_length'] < 3: reasons.append('fewer_than_three_pairs')
                if g['geometry_status'] != 'unique_at_numeric_tolerance': reasons.append('nonunique_rotation')
                result.update(numerical=numeric, geometry=g, numerical_usable=not reasons, numerical_exclusion_reasons=reasons)
                numeric_counts[mask + ':' + numeric['rmsd_status']] += 1; geometry_counts[mask + ':' + g['geometry_status']] += 1
                maximum = max(maximum, numeric['rmsd_rounding_error']); checked += 1
            else: result['numerical_exclusion_reasons'] = [record['status']]
            native_counts[mask + ':' + record['status']] += 1
            usable_counts[mask + ':' + ('usable' if result['numerical_usable'] else ';'.join(result['numerical_exclusion_reasons']))] += 1
            handle.write(json.dumps(result, separators=(',', ':')) + '\n')
            if index % 10000 == 0: print('Full new background numeric/geometry dispositions', index, '/', len(checkpoints), flush=True)
    assert dict(native_counts) == receipt['counts']
    verify(bindings)
    result = dict(status='complete_expanded_background_numeric_geometry_pending_independent_readback', plan_sha256=sha(plan_path), native_receipt_sha256=sha(root / 'receipt.json'),
                  directed_dispositions=len(checkpoints), numerically_checked_alignments=checked, counts=dict(native_counts), rmsd_status_counts=dict(numeric_counts), geometry_counts=dict(geometry_counts),
                  numerical_eligibility_counts=dict(usable_counts), maximum_rmsd_rounding_error=maximum, full_background_pairs=receipt['full_background_pairs'],
                  existing_catalog_pairs_pending_reuse=receipt['existing_catalog_pairs_pending_reuse'], source_hashes=bindings,
                  artifacts={name: sha(out / name) for name in ['disposition_measurements.jsonl.gz', ledger]}, scientific_eligibility=False, scope=plan['scope'])
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); run(parser.parse_args().plan)
