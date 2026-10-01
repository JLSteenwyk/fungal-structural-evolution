#!/usr/bin/env python3
"""Audit every new background state using separate mappings and quaternion fits."""
import argparse
import gzip
import json
import math
import re
from collections import Counter
from functools import lru_cache
from pathlib import Path
import numpy as np
from duplication_alignment_numeric_readback import load_pdb
from expanded_background_measurement_sources import load_native, inspect_native
from readback_background_alignment_geometry import verify_row
from readback_domain_triad_common_fits import quaternion_fit
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def numeric_fields(record, coords):
    """Separate raw-text, vectorized residue mapping and quaternion reconstruction."""
    lines = record['stdout'].splitlines(); m = record['metrics']
    summary = [line for line in lines if line.startswith('Aligned length=')]; assert len(summary) == 1
    match = re.fullmatch(r'Aligned length=\s*(\d+), RMSD=\s*([0-9.]+), Seq_ID=n_identical/n_aligned=\s*([0-9.]+)\s*', summary[0])
    assert match is not None
    n, native_rmsd, native_identity = int(match[1]), float(match[2]), float(match[3])
    markers = [i for i, line in enumerate(lines) if line.startswith('(":" denotes residue pairs')]; assert len(markers) == 1
    block = lines[markers[0] + 1:markers[0] + 4]; assert len(block) == 3
    strings = [block[0], block[2]]; assert block == [m['alignment_left'], m['alignment_marks'], m['alignment_right']]
    assert len(strings[0]) == len(strings[1]) == len(block[1])
    letters = [np.asarray(list(text)) for text in strings]; present = [a != '-' for a in letters]; paired = present[0] & present[1]
    assert int(paired.sum()) == n == m['aligned_length'] and n >= 1
    marks = np.asarray(list(block[1])); assert np.all(np.isin(marks[paired], [':', '.'])) and np.all(marks[~paired] == ' ')
    assert np.all(present[0] | present[1])
    indexes = [(np.cumsum(p) - 1)[paired] for p in present]
    lengths, scores = [], []
    for index, (text, c) in enumerate(zip(strings, coords), 1):
        assert text.replace('-', '') == c[0]
        l = [line for line in lines if line.startswith(f'Length of Structure_{index}:')]; assert len(l) == 1
        length = int(l[0].split(':', 1)[1].split()[0]); assert length == len(c[0]) == m['length_' + ('left' if index == 1 else 'right')]; lengths.append(length)
        s = [line for line in lines if line.startswith('TM-score=') and f'normalized by length of Structure_{index}:' in line]; assert len(s) == 1
        scores.append(float(s[0].split()[1]))
    x, y = [c[1][ix] for c, ix in zip(coords, indexes)]
    calculated, _, _ = quaternion_fit(x, y); error = abs(calculated - native_rmsd)
    identity = float(np.count_nonzero(letters[0][paired] == letters[1][paired]) / n)
    assert abs(identity - native_identity) <= .000501
    assert native_rmsd == m['rmsd'] and native_identity == m['sequence_identity'] and scores == [m['tm_left'], m['tm_right']]
    assert all(math.isfinite(v) for v in [native_rmsd, native_identity, *scores]) and native_rmsd >= 0 and all(0 <= v <= 1 for v in [native_identity, *scores])
    high = int(np.count_nonzero((coords[0][2][indexes[0]] >= 70) & (coords[1][2][indexes[1]] >= 70)))
    numeric = dict(rmsd_status='outside_printed_rounding' if error > .00501 else 'within_printed_rounding', aligned_length=n,
                   rmsd_recomputed=calculated, rmsd_native=native_rmsd, rmsd_rounding_error=error, sequence_identity_exact=identity,
                   tm_left_native=scores[0], tm_right_native=scores[1], coverage_left=n / lengths[0], coverage_right=n / lengths[1],
                   joint_plddt70_pairs=high, joint_plddt70_fraction=high / n)
    return numeric, x, y


def run(plan_path):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text())
    source, root, native, inputs, pairs, checkpoints, bindings, bundle = load_native(plan)
    bind(bindings, plan_path); native_plan_hash = sha(plan['native_plan'])
    ap = Path(plan['assessment_plan']); assessment_plan = json.loads(ap.read_text()); assert assessment_plan['native_plan'] == plan['native_plan']
    out = Path(assessment_plan['output']); rp = out / 'receipt.json'; r = json.loads(rp.read_text()); rh = sha(rp)
    assert r['status'] == 'complete_expanded_background_numeric_geometry_pending_independent_readback' and r['plan_sha256'] == sha(ap)
    assert r['native_receipt_sha256'] == sha(root / 'receipt.json') and r['scientific_eligibility'] is False
    expected_bindings = dict(bindings); expected_bindings.pop(str(plan_path))
    # Producer and reader pin the same scientific sources; each keeps its own plan.
    for path in plan.get('reader_only_pins', []): expected_bindings.pop(path)
    for path, digest in assessment_plan['pins'].items(): bind(expected_bindings, path, digest)
    bind(expected_bindings, ap); assert r['source_hashes'] == expected_bindings
    for path, digest in r['source_hashes'].items(): bind(bindings, path, digest)
    bind(bindings, rp, rh)
    for name, digest in r['artifacts'].items(): bind(bindings, out / name, digest)
    assert (out / 'full_background_work_partition.tsv').read_bytes() == (root / 'full_background_work_partition.tsv').read_bytes()
    seen = set(); native_counts, numeric_counts, geometry_counts, usable_counts = [Counter() for _ in range(4)]
    maximum_difference = maximum_curvature = maximum_error = 0.; checked = near = 0
    @lru_cache(maxsize=256)
    def coordinates(key): return load_pdb(inputs[key])
    with gzip.open(out / 'disposition_measurements.jsonl.gz', 'rt') as handle:
        for line in handle:
            row = json.loads(line); key = row['pair_key'], row['mask'], row['order']
            assert key in checkpoints and key not in seen; seen.add(key); pair, mask, order = key
            ends = pairs[pair] if order == 0 else pairs[pair][::-1]; rows = [inputs[(*end, mask)] for end in ends]; item = checkpoints[key]
            record = inspect_native(source, root, item, key, ends, rows, bundle); assert record['plan_sha256'] == native_plan_hash
            assert set(row) == {'pair_key', 'mask', 'order', 'model_left', 'version_left', 'model_right', 'version_right', 'checkpoint_path', 'checkpoint_sha256', 'native_status', 'native_elapsed_seconds', 'numerical', 'geometry', 'numerical_usable', 'numerical_exclusion_reasons'}
            assert [row['model_left'], row['version_left'], row['model_right'], row['version_right']] == [*ends[0], *ends[1]]
            assert row['checkpoint_path'] == str(root / item['path']) and row['checkpoint_sha256'] == item['sha256']
            assert row['native_status'] == record['status'] and row['native_elapsed_seconds'] == record['elapsed_seconds']
            if record['status'] == 'aligned':
                coords = [coordinates((*end, mask)) for end in ends]; numeric, x, y = numeric_fields(record, coords)
                assert set(row['numerical']) == set(numeric)
                for name, value in numeric.items():
                    actual = row['numerical'][name]
                    if isinstance(value, float): assert math.isfinite(actual) and math.isclose(actual, value, rel_tol=1e-9, abs_tol=1e-9), name
                    else: assert actual == value, name
                maximum_difference = max(maximum_difference, abs(row['numerical']['rmsd_recomputed'] - numeric['rmsd_recomputed']))
                err, boundary = verify_row(row['geometry'], x, y); maximum_curvature = max(maximum_curvature, err); near += boundary
                assert set(row['geometry']) == {'aligned_length', 'rank_left', 'rank_right', 'rms_width1_left', 'rms_width2_left', 'rms_width3_left', 'rms_width1_right', 'rms_width2_right', 'rms_width3_right', 'width2_to_width1_left', 'width2_to_width1_right', 'cross_s1', 'cross_s2', 'cross_s3', 'determinant_correction', 'minimum_rotation_curvature', 'relative_rotation_curvature', 'relative_numeric_tolerance', 'geometry_status'}
                reasons = []
                if numeric['rmsd_status'] != 'within_printed_rounding': reasons.append('rmsd_discrepancy')
                if numeric['aligned_length'] < 3: reasons.append('fewer_than_three_pairs')
                if row['geometry']['geometry_status'] != 'unique_at_numeric_tolerance': reasons.append('nonunique_rotation')
                numeric_counts[mask + ':' + numeric['rmsd_status']] += 1; geometry_counts[mask + ':' + row['geometry']['geometry_status']] += 1
                maximum_error = max(maximum_error, row['numerical']['rmsd_rounding_error']); checked += 1
                if mask == 'plddt70': assert numeric['joint_plddt70_pairs'] == numeric['aligned_length']
            else:
                assert row['numerical'] is None and row['geometry'] is None; reasons = [record['status']]
            usable = record['status'] == 'aligned' and not reasons
            assert row['numerical_exclusion_reasons'] == reasons and row['numerical_usable'] is usable
            native_counts[mask + ':' + record['status']] += 1; usable_counts[mask + ':' + ('usable' if usable else ';'.join(reasons))] += 1
            if len(seen) % 10000 == 0: print('Independent full new-background measurement dispositions', len(seen), '/', len(checkpoints), flush=True)
    summary = dict(directed_dispositions=len(seen), numerically_checked_alignments=checked, counts=dict(native_counts), rmsd_status_counts=dict(numeric_counts), geometry_counts=dict(geometry_counts),
                   numerical_eligibility_counts=dict(usable_counts), maximum_rmsd_rounding_error=maximum_error, full_background_pairs=native['full_background_pairs'], existing_catalog_pairs_pending_reuse=native['existing_catalog_pairs_pending_reuse'])
    assert seen == set(checkpoints) and dict(native_counts) == native['counts'] and all(r[k] == v for k, v in summary.items())
    verify(bindings)
    result = dict(status='passed_full_expanded_background_measurement_quaternion_readback', plan_sha256=sha(plan_path), producer_receipt_sha256=rh, **summary,
                  maximum_absolute_rmsd_difference=maximum_difference, maximum_scaled_quaternion_curvature_error=maximum_curvature, near_zero_quaternion_gaps=near,
                  source_hashes=bindings, scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['output']).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); run(parser.parse_args().plan)
