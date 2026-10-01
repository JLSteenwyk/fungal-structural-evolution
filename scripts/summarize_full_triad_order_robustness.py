#!/usr/bin/env python3
"""Retain full order grids and summarize strict all-order geometric qualification."""
import argparse
import csv
import gzip
import itertools
import json
import math
import shutil
from collections import Counter
from pathlib import Path
from full_triad_robustness_sources import load_sources, DEFINITIONS, NUMERIC_FIELDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def summarize(rows, screens):
    statuses = [r['fit_status'] for r in rows]; unique = all(s == 'computed_unique_at_numeric_tolerance' for s in statuses)
    ranges = {}
    for field in NUMERIC_FIELDS:
        numbers = [float(row[field]) for row in rows if row[field] != '']; assert all(math.isfinite(v) for v in numbers)
        ranges[field] = dict(available_orders=len(numbers), minimum=min(numbers) if numbers else None, maximum=max(numbers) if numbers else None,
                             span=max(numbers) - min(numbers) if numbers else None, mean=math.fsum(numbers) / len(numbers) if numbers else None)
    contrast = ranges['rmsd_ar_minus_br']
    direction = ('uncomputed_or_excluded_order' if contrast['available_orders'] != 8 else 'nonunique_order_fit' if not unique else
                 'positive' if contrast['minimum'] > 0 else 'negative' if contrast['maximum'] < 0 else 'includes_zero')
    result = dict(order_fit_statuses=statuses, all_orders_unique=unique, computed_orders=sum(s.startswith('computed_') for s in statuses),
                  order_triplet_hashes=[r['triples_sha256'] for r in rows], distinct_order_triplet_sets=len({r['triples_sha256'] for r in rows}),
                  order_source_exclusions=[r['source_exclusions'].split(';') if r['source_exclusions'] else [] for r in rows],
                  metric_ranges=ranges, strict_all_orders_contrast_direction=direction, screens={})
    for screen in screens:
        sid = screen['id']; flags = [int(r[sid + '_pass']) for r in rows]; core = [int(r[sid + '_core_pass']) for r in rows]; inherited = [int(r[sid + '_three_pair_pass']) for r in rows]
        assert all(v in [0, 1] for v in flags + core + inherited) and all(a == (b and c) for a, b, c in zip(flags, core, inherited))
        bits = sum(flag << i for i, flag in enumerate(flags))
        result['screens'][sid] = dict(order_pass_bits=bits, passing_orders=sum(flags), all_orders_pass=bits == 255, any_order_pass=bits != 0,
                                     core_order_pass_bits=sum(flag << i for i, flag in enumerate(core)), inherited_order_pass_bits=sum(flag << i for i, flag in enumerate(inherited)),
                                     order_exclusions=[r[sid + '_exclusions'].split(';') if r[sid + '_exclusions'] else [] for r in rows])
        if bits == 255: assert unique
    return result


def run(plan_path):
    plan = json.loads(Path(plan_path).read_text()); assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2 ** 30
    source, triads, contexts, work, bindings = load_sources(plan, plan_path); out = Path(plan['output']); out.mkdir(exist_ok=False, parents=True)
    counts = Counter(); passes = Counter(); any_passes = Counter(); maximum = {field: 0. for field in NUMERIC_FIELDS}; groups = rows = 0
    with gzip.open(source, 'rt') as f, gzip.open(out / 'triad_order_robustness.jsonl.gz', 'xt', compresslevel=1) as dest:
        reader = csv.DictReader(f, delimiter='\t')
        for number, triad in enumerate(triads, 1):
            for mask in ['full', 'plddt70']:
                by_definition = {d: [] for d in DEFINITIONS}
                for order in itertools.product([0, 1], repeat=3):
                    for definition in DEFINITIONS:
                        row = next(reader, None); assert row is not None
                        assert (row['triad_id'], row['mask'], row['mapping_definition']) == (triad['triad_id'], mask, definition)
                        assert tuple(int(row['order_' + edge]) for edge in ['ab', 'ar', 'br']) == order
                        assert [[row['model_' + role], int(row['version_' + role])] for role in ['a', 'b', 'reference']] == triad['models']
                        by_definition[definition].append(row); rows += 1
                for definition, states in by_definition.items():
                    assert len(states) == 8; value = summarize(states, plan['screens'])
                    record = dict(triad_id=triad['triad_id'], role_order=['a', 'b', 'reference'], models=triad['models'], mask=mask, mapping_definition=definition,
                                  order_bits_layout=[list(o) for o in itertools.product([0, 1], repeat=3)], **value)
                    dest.write(json.dumps(record, separators=(',', ':')) + '\n'); groups += 1
                    key = mask + ':' + definition; counts[key + ':' + value['strict_all_orders_contrast_direction']] += 1
                    counts[key + ':all_orders_unique'] += value['all_orders_unique']
                    for sid, screen in value['screens'].items(): passes[key + ':' + sid] += screen['all_orders_pass']; any_passes[key + ':' + sid] += screen['any_order_pass']
                    for field, stats in value['metric_ranges'].items():
                        if stats['span'] is not None: maximum[field] = max(maximum[field], stats['span'])
            if number % 1000 == 0: print('Full triad all-order groups', number, '/', len(triads), flush=True)
        assert next(reader, None) is None
    assert rows == plan['expected']['fit_rows'] and groups == 4 * len(triads) == plan['expected']['robustness_groups']
    verify(bindings)
    result = dict(status='complete_full_triad_all_order_robustness_pending_independent_readback', plan_sha256=sha(plan_path), correspondence_work_triads=len(triads),
                  fit_rows=rows, robustness_groups=groups, counts=dict(counts), all_order_screen_pass_counts=dict(passes), any_order_screen_pass_counts=dict(any_passes), maximum_metric_spans=maximum,
                  source_hashes=bindings, artifacts={'triad_order_robustness.jsonl.gz': sha(out / 'triad_order_robustness.jsonl.gz')}, scientific_eligibility=False, scope=plan['scope'])
    with (out / 'receipt.json').open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); run(p.parse_args().plan)


if __name__ == '__main__': main()
