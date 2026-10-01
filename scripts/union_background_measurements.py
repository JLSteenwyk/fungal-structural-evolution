#!/usr/bin/env python3
"""Merge every closed new and reused background state without changing exclusions."""
import argparse
import gzip
import json
import math
import shutil
from collections import Counter
from pathlib import Path
from background_measurement_union_sources import load_sources, NUMERIC_INTS, NUMERIC_FLOATS, GEOMETRY_INTS, GEOMETRY_FLOATS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text()); source, bindings = load_sources(plan, plan_path)
    if 'resources' in plan:
        assert shutil.disk_usage(Path(plan['output']).parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    new = {}
    with gzip.open(source['new_measurements'], 'rt') as handle:
        for line in handle:
            row = json.loads(line); key = row['pair_key'], row['mask'], row['order']; assert key not in new; new[key] = row
    assert set(new) == set(source['checkpoints'])
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False); ledger = out / source['ledger'].name; ledger.write_bytes(source['ledger'].read_bytes())
    seen = set(); owners, statuses, counts = Counter(), Counter(), Counter()
    with gzip.open(source['reuse_ledger'], 'rt') as old, gzip.open(out / 'background_measurement_dispositions.jsonl.gz', 'wt') as target:
        for line in old:
            row = json.loads(line); key = row['pair_key'], row['mask'], row['order']; pair, mask, order = key
            assert key not in seen and pair in source['full'] and mask in ['full', 'plddt70'] and order in [0, 1]
            ends = source['full'][pair]; assert ends == [(row['model_a'], row['version_a']), (row['model_b'], row['version_b'])]
            desired = ends if order == 0 else ends[::-1]
            if pair in source['new_pairs']:
                assert row['reuse_status'] == 'new_native_measurement_pending' and row['selected_source'] == ''
                m = new[key]; proof = source['checkpoints'][key]; path, digest, native_order = proof['path'], proof['sha256'], order
                owner = 'background_new_native'; numeric, geometry = m['numerical'], m['geometry']
                assert m['checkpoint_path'] == path and m['checkpoint_sha256'] == digest and m['native_status'] == proof['status']
                assert [m['model_left'], m['version_left'], m['model_right'], m['version_right']] == [*desired[0], *desired[1]]
            else:
                assert row['reuse_status'] == 'verified_identical_input_checkpoint_and_retained_disposition' and row['selected_source'] in ['reference_old', 'background_old']
                path, digest, native_order = row['source_checkpoint'], row['source_checkpoint_sha256'], row['source_order']; owner = row['selected_source']
                numeric, geometry = row['source_numeric'], row['source_geometry']
            assert native_order in [0, 1] and sha(path) == digest; bind(bindings, path, digest); native = json.loads(Path(path).read_text())
            assert (native['pair_key'], native['mask'], native['order']) == (pair, mask, native_order) and [(v['model_id'], v['version']) for v in native['inputs']] == desired
            status = native['status']; assert status in ['aligned', 'input_unavailable', 'native_error', 'parse_error', 'timeout']
            if owner == 'background_new_native':
                assert native['plan_sha256'] == source['native_plan_sha256'] and native['input_manifest_sha256'] == source['native_bundle'] and status == m['native_status']
            else: assert status == row['source_native_status']
            if status == 'aligned':
                assert native['returncode'] == 0 and numeric is not None and geometry is not None
                if owner != 'background_new_native':
                    assert (numeric['pair_key'], numeric['mask'], int(numeric['order'])) == (pair, mask, native_order)
                    assert (geometry['pair_key'], geometry['mask'], int(geometry['order'])) == (pair, mask, native_order) and geometry['rmsd_status'] == numeric['rmsd_status']
                n = {'rmsd_status': numeric['rmsd_status'], **{name: int(numeric[name]) for name in NUMERIC_INTS}, **{name: float(numeric[name]) for name in NUMERIC_FLOATS}}
                g = {'geometry_status': geometry['geometry_status'], **{name: int(geometry[name]) for name in GEOMETRY_INTS}, **{name: float(geometry[name]) for name in GEOMETRY_FLOATS}}
                assert n['rmsd_status'] in ['within_printed_rounding', 'outside_printed_rounding'] and g['geometry_status'] in ['unique_at_numeric_tolerance', 'degenerate_at_numeric_tolerance']
                assert all(math.isfinite(n[name]) for name in NUMERIC_FLOATS) and all(math.isfinite(g[name]) for name in GEOMETRY_FLOATS)
                assert n['rmsd_native'] == native['metrics']['rmsd'] and n['tm_left_native'] == native['metrics']['tm_left'] and n['tm_right_native'] == native['metrics']['tm_right']
                assert n['aligned_length'] == g['aligned_length'] == native['metrics']['aligned_length']; why = []
                if n['rmsd_status'] != 'within_printed_rounding': why.append('rmsd_discrepancy')
                if n['aligned_length'] < 3: why.append('fewer_than_three_pairs')
                if g['geometry_status'] != 'unique_at_numeric_tolerance': why.append('nonunique_rotation')
            else:
                assert numeric is None and geometry is None and 'metrics' not in native; n = g = None; why = [status]
            provenance = m if owner == 'background_new_native' else row
            assert provenance['numerical_exclusion_reasons'] == why and provenance['numerical_usable'] is (not why)
            exported = dict(pair_key=pair, model_a=ends[0][0], version_a=ends[0][1], model_b=ends[1][0], version_b=ends[1][1], mask=mask, order=order,
                directed_endpoints=[dict(model_id=end[0], version=end[1]) for end in desired], selected_source=owner, source_checkpoint=path, source_checkpoint_sha256=digest, source_order=native_order,
                source_plan_sha256=native['plan_sha256'], source_input_manifest_sha256=native['input_manifest_sha256'], source_native_status=status, source_numeric=numeric, source_geometry=geometry,
                numerical=n, geometry=g, native_metrics=native.get('metrics'), native_elapsed_seconds=native['elapsed_seconds'], native_returncode=native.get('returncode'),
                native_timeout_seconds=native.get('timeout_seconds'), native_error=native.get('error'), numerical_exclusion_reasons=why, numerical_usable=not why)
            target.write(json.dumps(exported, separators=(',', ':')) + '\n'); seen.add(key); owners[owner] += 1; statuses[mask + ':' + status] += 1; counts['usable' if not why else ';'.join(why)] += 1
            if len(seen) % 10000 == 0: print('Full background measurement union states', len(seen), '/', 4 * plan['full_pairs'], flush=True)
    assert seen == {(pair, mask, order) for pair in source['full'] for mask in ['full', 'plddt70'] for order in [0, 1]} and owners['background_new_native'] == 4 * plan['new_pairs']
    verify(bindings)
    result = dict(status='complete_full_background_measurement_union_pending_independent_readback', plan_sha256=sha(plan_path), full_pairs=plan['full_pairs'], directed_dispositions=len(seen),
                  source_dispositions=dict(owners), native_status_counts=dict(statuses), numerical_counts=dict(counts), source_hashes=bindings,
                  artifacts={name: sha(out / name) for name in ['background_measurement_dispositions.jsonl.gz', ledger.name]}, scientific_eligibility=False, scope=plan['scope'])
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); run(parser.parse_args().plan)
