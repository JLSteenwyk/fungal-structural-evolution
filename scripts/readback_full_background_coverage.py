#!/usr/bin/env python3
"""Independently reconstruct all background screens with SQL and decimal ceilings."""
import argparse
import csv
import gzip
import json
import sqlite3
import tempfile
from collections import Counter
from decimal import Decimal, ROUND_CEILING
from pathlib import Path
from full_background_coverage_sources import load, fields, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path, output):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text()); source, bindings = load(plan, plan_path); original_bindings = dict(bindings)
    out = Path(plan['output']); rp = out / 'receipt.json'; receipt = json.loads(rp.read_text()); rh = sha(rp)
    assert receipt['status'] == 'complete_full_background_coverage_pending_independent_readback' and receipt['plan_sha256'] == sha(plan_path) and receipt['scientific_eligibility'] is False
    assert receipt['source_hashes'] == original_bindings; bind(bindings, rp, rh)
    for name, digest in receipt['artifacts'].items(): bind(bindings, out / name, digest)
    verify(bindings)
    with tempfile.TemporaryDirectory(prefix='background-coverage-sql-', dir=out) as directory:
        db = sqlite3.connect(Path(directory) / 'original.sqlite'); db.execute('CREATE TABLE originals(pair TEXT,mask TEXT,ord INTEGER,record TEXT,PRIMARY KEY(pair,mask,ord))')
        states = 0
        with gzip.open(source['states'], 'rt') as handle:
            for line in handle:
                raw = json.loads(line); pair, mask, order = raw['pair_key'], raw['mask'], raw['order']
                assert pair in source['pairs'] and mask in ['full', 'plddt70'] and order in [0, 1]
                db.execute('INSERT INTO originals VALUES(?,?,?,?)', (pair, mask, order, line)); states += 1
        db.commit(); assert states == 4 * len(source['pairs'])
        flags, counts, exclusions = {}, Counter(), Counter()
        with (out / 'pair_mask_coverage.tsv').open() as handle:
            reader = csv.DictReader(handle, delimiter='\t'); assert reader.fieldnames == fields(plan['screens'])
            for actual in reader:
                pair, mask = actual['pair_key'], actual['mask']; key = pair, mask; assert key not in flags and pair in source['pairs'] and mask in ['full', 'plddt70']
                models = [source['models'][end] for end in source['pairs'][pair]]
                expected = dict(pair_key=pair, mask=mask, model_a=models[0]['model_id'], version_a=str(models[0]['version']), model_b=models[1]['model_id'], version_b=str(models[1]['version']))
                for side, model in zip(['a', 'b'], models):
                    for name in ['length', 'mean_ca_plddt', 'fraction_ca_plddt_below50']: expected[name + '_' + side] = str(model[name])
                eligible = []; not_usable = []
                for order in [0, 1]:
                    lookup = db.execute('SELECT record FROM originals WHERE pair=? AND mask=? AND ord=?', (pair, mask, order)).fetchone(); assert lookup is not None; raw = json.loads(lookup[0])
                    assert [(v['model_id'], v['version']) for v in raw['directed_endpoints']] == (source['pairs'][pair] if order == 0 else source['pairs'][pair][::-1])
                    n = raw['numerical']['aligned_length'] if raw['numerical'] is not None else None
                    assert (n is not None) == (raw['source_native_status'] == 'aligned') and raw['numerical_usable'] is (len(raw['numerical_exclusion_reasons']) == 0)
                    assert n is None or isinstance(n, int) and 0 < n <= min(v['length'] for v in models)
                    data = dict(selected_source=raw['selected_source'], source_order=raw['source_order'], source_checkpoint_sha256=raw['source_checkpoint_sha256'], native_status=raw['source_native_status'], numerical_usable=int(raw['numerical_usable']), numerical_exclusion_reasons=';'.join(raw['numerical_exclusion_reasons']), aligned_length=n if n is not None else '', original_coverage_a=n / models[0]['length'] if n is not None else '', original_coverage_b=n / models[1]['length'] if n is not None else '')
                    expected.update({f'order{order}_{name}': str(value) for name, value in data.items()})
                    if raw['numerical_usable']: eligible.append(n)
                    else: not_usable.append(f'order{order}_not_numerically_usable')
                flags[key] = {}
                for spec in plan['screens']:
                    minimum = max(int((Decimal(model['length']) * Decimal(str(spec['minimum_original_coverage']))).to_integral_value(rounding=ROUND_CEILING)) for model in models)
                    why = list(not_usable)
                    if any(n < spec['minimum_aligned_residues'] for n in eligible): why.append('short_alignment')
                    if any(n < minimum for n in eligible): why.append('low_original_protein_coverage')
                    name = spec['id']; expected[name + '_pass'] = str(int(not why)); expected[name + '_exclusions'] = ';'.join(why); flags[key][name] = not why
                    counts[mask + ':' + name] += not why
                    for reason in why: exclusions[mask + ':' + name + ':' + reason] += 1
                assert actual == expected, pair
                if len(flags) % 10000 == 0: print('Independent full background coverage rows', len(flags), '/', 2 * len(source['pairs']), flush=True)
        assert set(flags) == {(pair, mask) for pair in source['pairs'] for mask in ['full', 'plddt70']}
        both = {s['id']: sum(flags[p, 'full'][s['id']] and flags[p, 'plddt70'][s['id']] for p in source['pairs']) for s in plan['screens']}
        target = source['target']; summary = dict(source_inventory_models=len(source['models']), pair_models=source['pair_models'], pairs=len(source['pairs']), pair_mask_rows=len(flags), directed_source_states=states, screen_decisions=len(flags) * len(plan['screens']), screens=plan['screens'], pass_counts=dict(counts), both_masks_pass_counts=both, exclusion_counts=dict(exclusions), target_pairs=target['pairs'], target_pass_counts=target['pair_pass_counts'], target_both_masks_pass_counts=target['pair_both_masks_pass_counts'])
        assert all(receipt[name] == summary[name] for name in SUMMARY_FIELDS)
        expected_compare = {}
        for label, total, per_mask, joint in [('duplication_target', target['pairs'], target['pair_pass_counts'], target['pair_both_masks_pass_counts']), ('background', len(source['pairs']), dict(counts), both)]:
            for mask in ['full', 'plddt70', 'both']:
                for spec in plan['screens']:
                    key = label, mask, spec['id']; expected_compare[key] = dict(comparison_class=label, mask=mask, screen=spec['id'], total_physical_pairs=str(total), passed_physical_pairs=str(joint[spec['id']] if mask == 'both' else per_mask[mask + ':' + spec['id']]))
        with (out / 'target_background_physical_counts.tsv').open() as handle:
            for actual in csv.DictReader(handle, delimiter='\t'):
                key = actual['comparison_class'], actual['mask'], actual['screen']; assert expected_compare.pop(key) == actual
        assert not expected_compare; db.close()
    verify(bindings)
    result = dict(status='passed_full_background_coverage_sql_decimal_readback', plan_sha256=sha(plan_path), producer_receipt_sha256=rh, **summary, source_hashes=bindings, scientific_eligibility=False, scope=plan['scope'])
    with Path(output).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2)); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); run(args.plan, args.output)
