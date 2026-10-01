#!/usr/bin/env python3
"""Exercise full source I/O, exact coverage boundaries and independently rejected exports.

Prior native/geometry/union/target/journal proofs are explicitly synthetic
software fixture contracts. This does not prove physical measurements or
production closure and is not a pilot. Actual closed-source byte bindings,
full-grid handling, thresholds, missingness, export identity and SQL/decimal
readback run unchanged.
"""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha

SCREENS = [dict(id=f'n{n}_c{c}', minimum_aligned_residues=n, minimum_original_coverage=c / 100) for n in [30, 50] for c in [50, 70, 90]]


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True); path.write_text(json.dumps(value, indent=2) + '\n')


def table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n'); writer.writeheader(); writer.writerows(rows)


def zipped(path, rows):
    with gzip.open(path, 'wt') as handle:
        for row in rows: handle.write(json.dumps(row) + '\n')


def invoke(command, success=True):
    result = subprocess.run(command, capture_output=True, text=True)
    if success and result.returncode: raise RuntimeError(result.stderr + result.stdout)
    if not success: assert result.returncode != 0, 'Corrupt export accepted'


def setup(root):
    union, target = root / 'union', root / 'target'; union.mkdir(); target.mkdir()
    models = [dict(model_id=name, version=1 if name != 'B' else 2, length=n, mean_ca_plddt=90., fraction_ca_plddt_below50=0.) for name, n in [('A', 100), ('B', 101), ('C', 60), ('D', 120), ('E', 90)]]
    mp = root / 'models.jsonl.gz'; zipped(mp, models); lookup = {r['model_id']: r for r in models}
    pairs, states = [], []
    for b in ['B', 'C', 'D']:
        ends = [('A', 1), (b, lookup[b]['version'])]; pair = hashlib.sha256(json.dumps(ends, separators=(',', ':')).encode()).hexdigest()
        pairs.append(dict(pair_key=pair, model_a='A', version_a=1, model_b=b, version_b=lookup[b]['version']))
        for mask in ['full', 'plddt70']:
            for order in [0, 1]:
                state, why = 'aligned', []
                if b == 'B': n = (71 if order == 0 else 70) if mask == 'full' else 70
                elif b == 'C':
                    n = 50 if mask == 'full' else 30
                    if mask == 'plddt70' and order == 1: state, why = 'native_error', ['native_error']
                else:
                    n = 50 if order == 0 else 100
                    if mask == 'full' and order == 0: why = ['rmsd_discrepancy', 'nonunique_rotation']
                    if mask == 'plddt70': state, why = ('input_unavailable', ['input_unavailable']) if order == 0 else ('timeout', ['timeout'])
                directed = ends if order == 0 else ends[::-1]
                states.append(dict(pair_key=pair, mask=mask, order=order, directed_endpoints=[dict(model_id=e[0], version=e[1]) for e in directed], source_native_status=state,
                                   selected_source='reference_old' if b == 'B' else 'background_new_native', source_order=1 - order if b == 'B' else order, source_checkpoint_sha256='1' * 64,
                                   numerical=dict(aligned_length=n) if state == 'aligned' else None, numerical_exclusion_reasons=why, numerical_usable=not why))
    pairs.sort(key=lambda r: r['pair_key']); states.sort(key=lambda r: (r['pair_key'], r['mask'] != 'full', r['order']))
    table(union / 'full_background_work_partition.tsv', pairs); zipped(union / 'background_measurement_dispositions.jsonl.gz', states)
    up = root / 'union-plan.json'; write(up, dict(output=str(union)))
    summary = dict(full_pairs=3, directed_dispositions=12, source_dispositions=dict(Counter(r['selected_source'] for r in states)), native_status_counts=dict(Counter(r['mask'] + ':' + r['source_native_status'] for r in states)), numerical_counts=dict(Counter('usable' if r['numerical_usable'] else ';'.join(r['numerical_exclusion_reasons']) for r in states)))
    rp, ap = union / 'receipt.json', union / 'readback.json'
    write(rp, dict(status='complete_full_background_measurement_union_pending_independent_readback', plan_sha256=sha(up), **summary, scientific_eligibility=False, artifacts={name: sha(union / name) for name in ['full_background_work_partition.tsv', 'background_measurement_dispositions.jsonl.gz']}))
    write(ap, dict(status='passed_full_background_measurement_union_sql_readback', plan_sha256=sha(up), producer_receipt_sha256=sha(rp), **summary, scientific_eligibility=False))
    bindings = {str(p): sha(p) for p in [up, mp, rp, ap, union / 'full_background_work_partition.tsv', union / 'background_measurement_dispositions.jsonl.gz']}
    archive, completion = union / 'archive.json', root / 'union-completed.json'
    write(archive, dict(status='complete_verified_full_background_measurement_union_archive', services=[dict(synthetic_fixture_only=True)] * 2, source_hashes=bindings, summary=summary))
    write(completion, dict(status='complete_verified_full_background_measurement_union', **summary, scientific_eligibility=False, exact_process_journals_checked=2, full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive), bound_source_hashes=len(bindings), source_plan=str(up), source_plan_sha256=sha(up), producer_receipt=str(rp), producer_receipt_sha256=sha(rp), independent_readback=str(ap), independent_readback_sha256=sha(ap)))
    tp, sp = root / 'target-plan.json', root / 'screens.json'; write(tp, dict(output=str(target), screens=SCREENS)); write(sp, dict(screens=SCREENS))
    t = dict(status='complete_full_expanded_duplication_coverage_pending_independent_readback', plan_sha256=sha(tp), pairs=1, pair_mask_rows=2, screens=SCREENS, pair_pass_counts={m + ':' + s['id']: 0 for m in ['full', 'plddt70'] for s in SCREENS}, pair_both_masks_pass_counts={s['id']: 0 for s in SCREENS}); write(target / 'receipt.json', t)
    tc = root / 'target-completed.json'; write(tc, dict(status='complete_verified_full_expanded_pair_event_and_taxon_coverage_screens', services=[dict(synthetic_fixture_only=True)] * 2, scientific_eligibility=False, summary={k: t[k] for k in ['pairs', 'pair_mask_rows', 'pair_pass_counts', 'pair_both_masks_pass_counts']}, source_hashes={str(target / 'receipt.json'): sha(target / 'receipt.json'), str(tp): sha(tp)}))
    pp = root / 'plan.json'; write(pp, dict(measurement_completion=str(completion), measurement_plan=str(up), models=str(mp), target_completion=str(tc), target_plan=str(tp), screens_plan=str(sp), screens=SCREENS, output=str(root / 'coverage'), expected=dict(pairs=3, models=5, pair_models=4, target_pairs=1), resources=dict(minimum_free_disk_gib=0), pins={}, scope=__doc__))
    return pp, pairs


def run():
    with tempfile.TemporaryDirectory(prefix='full-background-coverage-fixture-') as directory:
        root = Path(directory); plan, pairs = setup(root); out = root / 'coverage'
        producer = [sys.executable, 'scripts/screen_full_background_coverage.py', '--plan', str(plan)]
        reader = [sys.executable, 'scripts/readback_full_background_coverage.py', '--plan', str(plan), '--output']
        invoke(producer); invoke(reader + [str(root / 'passed.json')])
        path = out / 'pair_mask_coverage.tsv'; original = list(csv.DictReader(path.open(), delimiter='\t')); original_bytes = path.read_bytes()
        receipt = json.loads((out / 'receipt.json').read_text()); receipt_bytes = (out / 'receipt.json').read_bytes()
        by = {(r['model_b'], r['mask']): r for r in original}
        assert by['B', 'full']['n30_c70_pass'] == by['B', 'plddt70']['n30_c70_pass'] == '0'
        assert by['B', 'full']['order0_original_coverage_b'] == str(71 / 101) and by['B', 'full']['order1_original_coverage_b'] == str(70 / 101)
        assert by['B', 'plddt70']['n30_c70_exclusions'] == 'low_original_protein_coverage'
        assert by['C', 'full']['n50_c50_pass'] == '1' and by['C', 'plddt70']['n50_c50_exclusions'] == 'order1_not_numerically_usable;short_alignment;low_original_protein_coverage'
        assert by['D', 'full']['n30_c50_exclusions'] == 'order0_not_numerically_usable' and by['D', 'plddt70']['n30_c50_exclusions'] == 'order0_not_numerically_usable;order1_not_numerically_usable'
        assert receipt['pairs'] == 3 and receipt['pair_mask_rows'] == 6 and receipt['screen_decisions'] == 36
        assert receipt['source_inventory_models'] == 5 and receipt['pair_models'] == 4
        assert receipt['both_masks_pass_counts'] == {s['id']: int(s['minimum_original_coverage'] == 0.5) for s in SCREENS}
        rejected = []
        labels = ['missing_row', 'duplicate_row', 'wrong_original_length', 'wrong_masked_denominator', 'boundary_pass', 'cleared_numerical_flag', 'cleared_native_error', 'swapped_model_roles', 'changed_original_order', 'changed_source', 'changed_checkpoint', 'changed_summary', 'changed_target_comparison']
        comparison = out / 'target_background_physical_counts.tsv'; comparison_bytes = comparison.read_bytes()
        for label in labels:
            rows = copy.deepcopy(original); r = copy.deepcopy(receipt); choose = {(x['model_b'], x['mask']): x for x in rows}
            if label == 'missing_row': rows.pop()
            elif label == 'duplicate_row': rows.append(copy.deepcopy(rows[0]))
            elif label == 'wrong_original_length': choose['B', 'plddt70']['length_b'] = '70'
            elif label == 'wrong_masked_denominator': choose['B', 'plddt70']['order0_original_coverage_b'] = '1.0'
            elif label == 'boundary_pass': choose['B', 'full'].update(n30_c70_pass='1', n30_c70_exclusions='')
            elif label == 'cleared_numerical_flag': choose['D', 'full']['order0_numerical_exclusion_reasons'] = ''
            elif label == 'cleared_native_error': choose['C', 'plddt70']['order1_native_status'] = 'aligned'
            elif label == 'swapped_model_roles': rows[0]['model_a'], rows[0]['model_b'] = rows[0]['model_b'], rows[0]['model_a']
            elif label == 'changed_original_order': rows[0]['order0_source_order'] = str(1 - int(rows[0]['order0_source_order']))
            elif label == 'changed_source': rows[0]['order0_selected_source'] = 'other_source'
            elif label == 'changed_checkpoint': rows[0]['order0_source_checkpoint_sha256'] = '0' * 64
            elif label == 'changed_summary': r['pass_counts']['full:n30_c50'] += 1
            else:
                items = list(csv.DictReader(comparison.open(), delimiter='\t')); items[0]['total_physical_pairs'] = '999'; table(comparison, items)
            table(path, rows); r['artifacts'][path.name] = sha(path); r['artifacts'][comparison.name] = sha(comparison); write(out / 'receipt.json', r)
            invoke(reader + [str(root / (label + '.json'))], success=False); rejected.append(label); comparison.write_bytes(comparison_bytes)
        path.write_bytes(original_bytes); (out / 'receipt.json').write_bytes(receipt_bytes)
        return dict(status='passed_full_background_coverage_software_cases', source_inventory_models=5, pair_models=4, broader_inventory_retained=True, physical_pairs=3, pair_mask_rows=6, directed_states=12, screen_decisions=36, exact_integer_fraction_vs_decimal_ceiling_boundary=True, confidence_mask_uses_original_lengths=True, both_orders_and_both_mask_intersections_checked=True, rejected_rehashed_exports=rejected, synthetic_prior_measurement_target_and_journal_contracts=True, production_native_or_journal_correctness_tested=False, scope=__doc__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path); args = parser.parse_args(); result = run(); result['script_sha256'] = sha(__file__)
    if args.output:
        with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
