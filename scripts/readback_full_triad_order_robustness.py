#!/usr/bin/env python3
"""Independently aggregate the complete fit/order grid using SQLite."""
import argparse
import csv
import gzip
import json
import math
import sqlite3
from collections import Counter
from pathlib import Path
from full_triad_robustness_sources import load_sources, DEFINITIONS, NUMERIC_FIELDS, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path, output):
    plan_path, output = Path(plan_path), Path(output); assert not output.exists()
    plan = json.loads(plan_path.read_text()); source, triads, _, _, bindings = load_sources(plan, plan_path)
    root = Path(plan['output']); rp = root / 'receipt.json'; receipt = json.loads(rp.read_text())
    assert receipt['status'] == 'complete_full_triad_all_order_robustness_pending_independent_readback'
    assert receipt['plan_sha256'] == sha(plan_path) and receipt['source_hashes'] == bindings
    assert receipt['scope'] == plan['scope'] and receipt['scientific_eligibility'] is False
    assert set(receipt['artifacts']) == {'triad_order_robustness.jsonl.gz'}
    bind(bindings, rp)
    for name, digest in receipt['artifacts'].items(): bind(bindings, root / name, digest)
    verify(bindings)
    database = root / 'independent_order_readback.sqlite'; assert not database.exists()
    db = sqlite3.connect(database); db.execute('PRAGMA journal_mode=OFF'); db.execute('PRAGMA synchronous=OFF')
    numeric_columns = ','.join('n' + str(i) + ' REAL' for i in range(len(NUMERIC_FIELDS)))
    screen_columns = ','.join('s' + str(i) + '_' + suffix + (' TEXT' if suffix == 'ex' else ' INT') for i in range(len(plan['screens'])) for suffix in ['pass', 'core', 'inherited', 'ex'])
    db.execute('CREATE TABLE fits(t TEXT,m TEXT,d TEXT,o INT,status TEXT,h TEXT,ex TEXT,' + numeric_columns + ',' + screen_columns + ',PRIMARY KEY(t,m,d,o))')
    catalog = {t['triad_id']: t for t in triads}; assert len(catalog) == len(triads)
    rows = 0; batch = []; slots = 7 + len(NUMERIC_FIELDS) + 4 * len(plan['screens'])
    with gzip.open(source, 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            triad = catalog[row['triad_id']]
            assert [[row['model_' + role], int(row['version_' + role])] for role in ['a', 'b', 'reference']] == triad['models']
            assert row['mask'] in ['full', 'plddt70'] and row['mapping_definition'] in DEFINITIONS
            orders = [int(row['order_' + edge]) for edge in ['ab', 'ar', 'br']]; assert all(x in [0, 1] for x in orders)
            order_index = 4 * orders[0] + 2 * orders[1] + orders[2]
            numbers = [None if row[field] == '' else float(row[field]) for field in NUMERIC_FIELDS]
            assert all(value is None or math.isfinite(value) for value in numbers)
            values = [row['triad_id'], row['mask'], row['mapping_definition'], order_index, row['fit_status'], row['triples_sha256'], row['source_exclusions'], *numbers]
            for screen in plan['screens']:
                sid = screen['id']; flags = [int(row[sid + suffix]) for suffix in ['_pass', '_core_pass', '_three_pair_pass']]
                assert all(x in [0, 1] for x in flags) and flags[0] == (flags[1] and flags[2])
                values += [*flags, row[sid + '_exclusions']]
            batch.append(values); rows += 1
            if len(batch) == 10000:
                db.executemany('INSERT INTO fits VALUES(' + ','.join(['?'] * slots) + ')', batch); batch.clear()
                print('Independent SQL fit ingestion', rows, '/', plan['expected']['fit_rows'], flush=True)
    if batch: db.executemany('INSERT INTO fits VALUES(' + ','.join(['?'] * slots) + ')', batch)
    db.commit(); assert rows == plan['expected']['fit_rows']
    assert db.execute('SELECT COUNT(*) FROM (SELECT t,m,d FROM fits GROUP BY t,m,d)').fetchone()[0] == plan['expected']['robustness_groups']
    assert db.execute('SELECT COUNT(*) FROM (SELECT t,m,d FROM fits GROUP BY t,m,d HAVING COUNT(*)!=8 OR MIN(o)!=0 OR MAX(o)!=7)').fetchone()[0] == 0
    numeric_sql = ','.join(f'COUNT(n{i}),MIN(n{i}),MAX(n{i}),MAX(n{i})-MIN(n{i}),AVG(n{i})' for i in range(len(NUMERIC_FIELDS)))
    counts, passes, any_passes = Counter(), Counter(), Counter(); maximum = {field: 0. for field in NUMERIC_FIELDS}; groups = 0; max_difference = 0.
    with gzip.open(root / 'triad_order_robustness.jsonl.gz', 'rt') as handle:
        for triad in triads:
            for mask in ['full', 'plddt70']:
                for definition in DEFINITIONS:
                    line = next(handle, None); assert line is not None
                    actual = json.loads(line); key = (triad['triad_id'], mask, definition)
                    assert (actual['triad_id'], actual['mask'], actual['mapping_definition']) == key
                    expected = dict(triad_id=key[0], role_order=['a', 'b', 'reference'], models=triad['models'], mask=mask, mapping_definition=definition,
                                    order_bits_layout=[[i // 4, (i // 2) % 2, i % 2] for i in range(8)])
                    raw = list(db.execute('SELECT status,h,ex FROM fits WHERE t=? AND m=? AND d=? ORDER BY o', key)); assert len(raw) == 8
                    statuses = [x[0] for x in raw]; unique = all(s == 'computed_unique_at_numeric_tolerance' for s in statuses)
                    expected.update(order_fit_statuses=statuses, all_orders_unique=unique, computed_orders=sum(s.startswith('computed_') for s in statuses),
                                    order_triplet_hashes=[x[1] for x in raw], distinct_order_triplet_sets=len({x[1] for x in raw}),
                                    order_source_exclusions=[x[2].split(';') if x[2] else [] for x in raw])
                    aggregate = db.execute('SELECT ' + numeric_sql + ' FROM fits WHERE t=? AND m=? AND d=?', key).fetchone()
                    ranges = {}
                    assert set(actual['metric_ranges']) == set(NUMERIC_FIELDS)
                    for i, field in enumerate(NUMERIC_FIELDS):
                        values = aggregate[5 * i:5 * (i + 1)]; stats = dict(zip(['available_orders', 'minimum', 'maximum', 'span', 'mean'], values))
                        observed = actual['metric_ranges'][field]; assert set(observed) == set(stats)
                        for label, value in stats.items():
                            if value is None: assert observed[label] is None
                            elif label == 'available_orders': assert type(observed[label]) is int and observed[label] == value
                            else:
                                assert type(observed[label]) in [int, float] and math.isfinite(observed[label])
                                difference = abs(observed[label] - value); max_difference = max(max_difference, difference)
                                assert difference <= plan['numeric_readback_absolute_tolerance'], (key, field, label, difference)
                        ranges[field] = observed
                        if values[3] is not None: maximum[field] = max(maximum[field], values[3])
                    expected['metric_ranges'] = ranges
                    contrast = aggregate[5 * NUMERIC_FIELDS.index('rmsd_ar_minus_br'):5 * (NUMERIC_FIELDS.index('rmsd_ar_minus_br') + 1)]
                    direction = 'uncomputed_or_excluded_order' if contrast[0] < 8 else 'nonunique_order_fit' if not unique else 'positive' if contrast[1] > 0 else 'negative' if contrast[2] < 0 else 'includes_zero'
                    expected['strict_all_orders_contrast_direction'] = direction; expected['screens'] = {}
                    group_key = mask + ':' + definition; counts[group_key + ':' + direction] += 1; counts[group_key + ':all_orders_unique'] += unique
                    for i, screen in enumerate(plan['screens']):
                        sql = f'SELECT SUM(s{i}_pass << o),SUM(s{i}_pass),MIN(s{i}_pass),MAX(s{i}_pass),SUM(s{i}_core << o),SUM(s{i}_inherited << o) FROM fits WHERE t=? AND m=? AND d=?'
                        bits, n, all_pass, any_pass, core, inherited = db.execute(sql, key).fetchone()
                        exclusions = [x[0].split(';') if x[0] else [] for x in db.execute(f'SELECT s{i}_ex FROM fits WHERE t=? AND m=? AND d=? ORDER BY o', key)]
                        sid = screen['id']; expected['screens'][sid] = dict(order_pass_bits=bits, passing_orders=n, all_orders_pass=bool(all_pass), any_order_pass=bool(any_pass),
                                                                         core_order_pass_bits=core, inherited_order_pass_bits=inherited, order_exclusions=exclusions)
                        if all_pass: assert unique
                        passes[group_key + ':' + sid] += all_pass; any_passes[group_key + ':' + sid] += any_pass
                    assert actual == expected and type(actual['all_orders_unique']) is bool
                    assert all(type(s['all_orders_pass']) is type(s['any_order_pass']) is bool for s in actual['screens'].values())
                    groups += 1
                    if groups % 10000 == 0: print('Independent SQL all-order groups', groups, '/', plan['expected']['robustness_groups'], flush=True)
        assert next(handle, None) is None
    summary = dict(correspondence_work_triads=len(triads), fit_rows=rows, robustness_groups=groups, counts=dict(counts), all_order_screen_pass_counts=dict(passes),
                   any_order_screen_pass_counts=dict(any_passes), maximum_metric_spans=maximum)
    assert set(summary) == set(SUMMARY_FIELDS) and groups == plan['expected']['robustness_groups']
    assert all(receipt[k] == v for k, v in summary.items()); verify(bindings)
    db.close(); database.unlink()
    result = dict(status='passed_full_triad_all_order_robustness_sql_readback', plan_sha256=sha(plan_path), producer_receipt_sha256=sha(rp), checker_sha256=sha(__file__),
                  **summary, maximum_absolute_numeric_difference=max_difference, source_hashes=bindings, scientific_eligibility=False,
                  scope='All865792 source fit dispositions independently indexed by ordered physical triad/mask/core/order; all108224 groups, every missing value, range/sign/order bitmap/status/hash/exclusion and every summary reconstructed with SQL. No producer summary algorithm imported; shares closed-source I/O only. Strict numeric sign is descriptive, not calibrated evolutionary asymmetry; original parent/orthology eligibility is a later full-context linkage.')
    with output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); run(args.plan, args.output)
