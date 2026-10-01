#!/usr/bin/env python3
"""Independently rebuild every full background union field from closed original sources."""
import argparse
import gzip
import json
import math
import sqlite3
import tempfile
from collections import Counter
from pathlib import Path
from background_measurement_union_sources import load_sources, NUMERIC_INTS, NUMERIC_FLOATS, GEOMETRY_INTS, GEOMETRY_FLOATS, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path, output):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text()); source, bindings = load_sources(plan, plan_path)
    out = Path(plan['output']); rp = out / 'receipt.json'; receipt = json.loads(rp.read_text()); rh = sha(rp)
    assert receipt['status'] == 'complete_full_background_measurement_union_pending_independent_readback' and receipt['plan_sha256'] == sha(plan_path) and receipt['scientific_eligibility'] is False
    initial_bindings = dict(bindings); bind(bindings, rp, rh)
    for name, digest in receipt['artifacts'].items(): bind(bindings, out / name, digest)
    verify(bindings); assert (out / source['ledger'].name).read_bytes() == source['ledger'].read_bytes()
    with tempfile.TemporaryDirectory(prefix='independent-background-union-', dir=out) as directory:
        db = sqlite3.connect(Path(directory) / 'source.sqlite'); db.execute('CREATE TABLE states(pair TEXT,mask TEXT,ord INTEGER,record TEXT,PRIMARY KEY(pair,mask,ord))')
        keys = set()
        with gzip.open(source['new_measurements'], 'rt') as handle:
            for line in handle:
                row = json.loads(line); key = row['pair_key'], row['mask'], row['order']; assert key in source['checkpoints'] and key not in keys; keys.add(key)
                db.execute('INSERT INTO states VALUES(?,?,?,?)', (*key, line))
        db.commit(); assert keys == set(source['checkpoints'])
        seen = set(); owners, statuses, counts = Counter(), Counter(), Counter()
        with gzip.open(source['reuse_ledger'], 'rt') as original, gzip.open(out / 'background_measurement_dispositions.jsonl.gz', 'rt') as exported:
            for line in original:
                old = json.loads(line); pair, mask, order = old['pair_key'], old['mask'], old['order']; key = pair, mask, order
                assert key not in seen and pair in source['full'] and mask in ['full', 'plddt70'] and order in [0, 1]
                models = source['full'][pair]; assert models == [(old['model_a'], old['version_a']), (old['model_b'], old['version_b'])]
                directed = models if order == 0 else list(reversed(models)); lookup = db.execute('SELECT record FROM states WHERE pair=? AND mask=? AND ord=?', key).fetchone()
                if lookup is not None:
                    assert pair in source['new_pairs'] and old['reuse_status'] == 'new_native_measurement_pending' and not old['selected_source']
                    measurement = json.loads(lookup[0]); proof = source['checkpoints'][key]
                    label, path, digest, actual_order = 'background_new_native', proof['path'], proof['sha256'], order
                    numeric, geometric = measurement['numerical'], measurement['geometry']; expected_flags = measurement
                    assert measurement['checkpoint_path'] == path and measurement['checkpoint_sha256'] == digest and measurement['native_status'] == proof['status']
                    assert [(measurement['model_left'], measurement['version_left']), (measurement['model_right'], measurement['version_right'])] == directed
                else:
                    assert pair not in source['new_pairs'] and old['reuse_status'] == 'verified_identical_input_checkpoint_and_retained_disposition'
                    label, path, digest, actual_order = [old[name] for name in ['selected_source', 'source_checkpoint', 'source_checkpoint_sha256', 'source_order']]
                    assert label in ['reference_old', 'background_old']; numeric, geometric = old['source_numeric'], old['source_geometry']; expected_flags = old
                assert actual_order in [0, 1] and sha(path) == digest; bind(bindings, path, digest); bind(initial_bindings, path, digest)
                raw = json.loads(Path(path).read_text()); state = raw['status']
                assert (raw['pair_key'], raw['mask'], raw['order']) == (pair, mask, actual_order)
                assert [(r['model_id'], r['version']) for r in raw['inputs']] == directed
                if lookup is not None: assert state == proof['status'] and raw['plan_sha256'] == source['native_plan_sha256'] and raw['input_manifest_sha256'] == source['native_bundle']
                else: assert state == old['source_native_status']
                if state == 'aligned':
                    assert raw['returncode'] == 0 and numeric is not None and geometric is not None
                    if lookup is None:
                        assert numeric['pair_key'] == geometric['pair_key'] == pair and numeric['mask'] == geometric['mask'] == mask
                        assert int(numeric['order']) == int(geometric['order']) == actual_order and geometric['rmsd_status'] == numeric['rmsd_status']
                    number = {name: int(numeric[name]) for name in NUMERIC_INTS}; number.update({name: float(numeric[name]) for name in NUMERIC_FLOATS}); number['rmsd_status'] = numeric['rmsd_status']
                    shape = {name: int(geometric[name]) for name in GEOMETRY_INTS}; shape.update({name: float(geometric[name]) for name in GEOMETRY_FLOATS}); shape['geometry_status'] = geometric['geometry_status']
                    assert all(math.isfinite(value) for name, value in number.items() if name != 'rmsd_status') and all(math.isfinite(value) for name, value in shape.items() if name != 'geometry_status')
                    assert number['aligned_length'] == shape['aligned_length'] == raw['metrics']['aligned_length']
                    assert number['rmsd_status'] in ['within_printed_rounding', 'outside_printed_rounding'] and shape['geometry_status'] in ['unique_at_numeric_tolerance', 'degenerate_at_numeric_tolerance']
                    conditions = [('rmsd_discrepancy', number['rmsd_status'] == 'outside_printed_rounding'), ('fewer_than_three_pairs', number['aligned_length'] < 3), ('nonunique_rotation', shape['geometry_status'] == 'degenerate_at_numeric_tolerance')]
                    reasons = [name for name, value in conditions if value]
                    assert number['rmsd_native'] == raw['metrics']['rmsd'] and number['tm_left_native'] == raw['metrics']['tm_left'] and number['tm_right_native'] == raw['metrics']['tm_right']
                else:
                    assert state in ['input_unavailable', 'native_error', 'parse_error', 'timeout'] and numeric is None and geometric is None and 'metrics' not in raw
                    number = shape = None; reasons = [state]
                assert reasons == expected_flags['numerical_exclusion_reasons'] and expected_flags['numerical_usable'] is (len(reasons) == 0)
                expected = dict(pair_key=pair, model_a=models[0][0], version_a=models[0][1], model_b=models[1][0], version_b=models[1][1], mask=mask, order=order,
                    directed_endpoints=[dict(model_id=r['model_id'], version=r['version']) for r in raw['inputs']], selected_source=label, source_checkpoint=path, source_checkpoint_sha256=digest, source_order=actual_order,
                    source_plan_sha256=raw['plan_sha256'], source_input_manifest_sha256=raw['input_manifest_sha256'], source_native_status=state, source_numeric=numeric, source_geometry=geometric,
                    numerical=number, geometry=shape, native_metrics=raw.get('metrics'), native_elapsed_seconds=raw['elapsed_seconds'], native_returncode=raw.get('returncode'),
                    native_timeout_seconds=raw.get('timeout_seconds'), native_error=raw.get('error'), numerical_exclusion_reasons=reasons, numerical_usable=len(reasons) == 0)
                actual = json.loads(next(exported)); assert actual == expected, pair
                seen.add(key); owners[label] += 1; statuses[mask + ':' + state] += 1; counts['usable' if not reasons else ';'.join(reasons)] += 1
                if len(seen) % 10000 == 0: print('Independent full background measurement union states', len(seen), '/', 4 * plan['full_pairs'], flush=True)
            assert next(exported, None) is None
        assert seen == {(pair, mask, order) for pair in source['full'] for mask in ['full', 'plddt70'] for order in [0, 1]}
        summary = dict(full_pairs=len(source['full']), directed_dispositions=len(seen), source_dispositions=dict(owners), native_status_counts=dict(statuses), numerical_counts=dict(counts))
        assert all(receipt[name] == summary[name] for name in SUMMARY_FIELDS) and owners['background_new_native'] == 4 * plan['new_pairs']
        assert receipt['source_hashes'] == initial_bindings; db.close()
    verify(bindings)
    result = dict(status='passed_full_background_measurement_union_sql_readback', plan_sha256=sha(plan_path), producer_receipt_sha256=rh, **summary, source_hashes=bindings, scientific_eligibility=False, scope=plan['scope'])
    with Path(output).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); run(args.plan, args.output)
