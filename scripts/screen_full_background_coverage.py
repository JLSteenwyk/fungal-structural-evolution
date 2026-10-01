#!/usr/bin/env python3
"""Screen every full background pair under the unchanged duplication coverage policy."""
import argparse
import csv
import gzip
import json
import shutil
from collections import Counter
from fractions import Fraction
from pathlib import Path
from full_background_coverage_sources import load, fields
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text()); source, bindings = load(plan, plan_path)
    out = Path(plan['output']); assert shutil.disk_usage(out.parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30; out.mkdir(exist_ok=False)
    counts, exclusions, flags = Counter(), Counter(), {}; states = 0
    with gzip.open(source['states'], 'rt') as native, (out / 'pair_mask_coverage.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fields(plan['screens']), delimiter='\t', lineterminator='\n'); writer.writeheader()
        for pair, ends in sorted(source['pairs'].items()):
            originals = [source['models'][end] for end in ends]
            for mask in ['full', 'plddt70']:
                row = dict(pair_key=pair, mask=mask, model_a=ends[0][0], version_a=ends[0][1], model_b=ends[1][0], version_b=ends[1][1])
                for side, model in zip(['a', 'b'], originals):
                    for name in ['length', 'mean_ca_plddt', 'fraction_ca_plddt_below50']: row[name + '_' + side] = model[name]
                usable_lengths = []; rejected_orders = []
                for order in [0, 1]:
                    raw = json.loads(next(native)); assert (raw['pair_key'], raw['mask'], raw['order']) == (pair, mask, order)
                    assert [(v['model_id'], v['version']) for v in raw['directed_endpoints']] == (ends if order == 0 else ends[::-1])
                    numerical = raw['numerical']; n = numerical['aligned_length'] if numerical is not None else None
                    assert (raw['source_native_status'] == 'aligned') == (n is not None)
                    assert raw['numerical_usable'] is (not raw['numerical_exclusion_reasons'])
                    assert n is None or isinstance(n, int) and 0 < n <= min(v['length'] for v in originals)
                    values = dict(selected_source=raw['selected_source'], source_order=raw['source_order'], source_checkpoint_sha256=raw['source_checkpoint_sha256'], native_status=raw['source_native_status'], numerical_usable=int(raw['numerical_usable']), numerical_exclusion_reasons=';'.join(raw['numerical_exclusion_reasons']), aligned_length=n if n is not None else '', original_coverage_a=n / originals[0]['length'] if n is not None else '', original_coverage_b=n / originals[1]['length'] if n is not None else '')
                    row.update({f'order{order}_{k}': v for k, v in values.items()}); states += 1
                    if raw['numerical_usable']: usable_lengths.append(n)
                    else: rejected_orders.append(f'order{order}_not_numerically_usable')
                flags[pair, mask] = {}
                for spec in plan['screens']:
                    why = list(rejected_orders); cutoff = Fraction(str(spec['minimum_original_coverage']))
                    if usable_lengths and min(usable_lengths) < spec['minimum_aligned_residues']: why.append('short_alignment')
                    if any(n * cutoff.denominator < original['length'] * cutoff.numerator for n in usable_lengths for original in originals): why.append('low_original_protein_coverage')
                    name = spec['id']; row[name + '_pass'] = int(not why); row[name + '_exclusions'] = ';'.join(why); flags[pair, mask][name] = not why
                    counts[mask + ':' + name] += not why
                    for reason in why: exclusions[mask + ':' + name + ':' + reason] += 1
                writer.writerow(row)
            if len(flags) % 10000 == 0: print('Full background pair/mask coverage', len(flags), '/', 2 * len(source['pairs']), flush=True)
        assert next(native, None) is None
    assert states == 4 * len(source['pairs']) and len(flags) == 2 * len(source['pairs'])
    both = {s['id']: sum(flags[p, 'full'][s['id']] and flags[p, 'plddt70'][s['id']] for p in source['pairs']) for s in plan['screens']}
    target = source['target']; compare = out / 'target_background_physical_counts.tsv'
    with compare.open('w') as handle:
        writer = csv.DictWriter(handle, ['comparison_class', 'mask', 'screen', 'total_physical_pairs', 'passed_physical_pairs'], delimiter='\t', lineterminator='\n'); writer.writeheader()
        for label, total, per_mask, joint in [('duplication_target', target['pairs'], target['pair_pass_counts'], target['pair_both_masks_pass_counts']), ('background', len(source['pairs']), dict(counts), both)]:
            for mask in ['full', 'plddt70', 'both']:
                for spec in plan['screens']:
                    writer.writerow(dict(comparison_class=label, mask=mask, screen=spec['id'], total_physical_pairs=total, passed_physical_pairs=joint[spec['id']] if mask == 'both' else per_mask[mask + ':' + spec['id']]))
    verify(bindings)
    result = dict(status='complete_full_background_coverage_pending_independent_readback', plan_sha256=sha(plan_path), source_inventory_models=len(source['models']), pair_models=source['pair_models'], pairs=len(source['pairs']), pair_mask_rows=len(flags), directed_source_states=states, screen_decisions=len(flags) * len(plan['screens']), screens=plan['screens'], pass_counts=dict(counts), both_masks_pass_counts=both, exclusion_counts=dict(exclusions), target_pairs=target['pairs'], target_pass_counts=target['pair_pass_counts'], target_both_masks_pass_counts=target['pair_both_masks_pass_counts'], source_hashes=bindings, artifacts={name: sha(out / name) for name in ['pair_mask_coverage.tsv', compare.name]}, scientific_eligibility=False, scope=plan['scope'])
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2)); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); run(parser.parse_args().plan)
