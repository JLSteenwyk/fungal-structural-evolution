#!/usr/bin/env python3
"""Exercise full order/coverage exports, exact boundaries and corrupted-data rejection."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')


def table(path, rows):
    with path.open('w') as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)


def run(command):
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode: raise RuntimeError(result.stdout + result.stderr)
    return result


def main():
    with tempfile.TemporaryDirectory(prefix='reference-order-coverage-fixture-') as temp:
        root = Path(temp); inventory = root / 'inventory'; union = root / 'union'
        inventory.mkdir(); union.mkdir()
        lengths = dict(A=100, B=200, C=101, D=60, E=101, F=100)
        versions = dict(A=6, B=10, C=6, D=10, E=10, F=6)
        models = [dict(model_id=m, version=versions[m], length=n, sequence_sha256=hashlib.sha256(('A' * n).encode()).hexdigest(), sha256=hashlib.sha256(m.encode()).hexdigest()) for m, n in lengths.items()]
        (inventory / 'models.jsonl').write_text(''.join(json.dumps(m) + '\n' for m in models))
        pairs = []
        for a, b in [('B', 'A'), ('A', 'C'), ('B', 'C'), ('A', 'D'), ('E', 'F')]:
            ends = [(a, versions[a]), (b, versions[b])]
            pairs.append(dict(pair_key=hashlib.sha256(json.dumps(sorted(ends), separators=(',', ':')).encode()).hexdigest(), model_a=a, version_a=versions[a], model_b=b, version_b=versions[b], work_disposition='additional_pair'))
        table(inventory / 'model_pairs.tsv', pairs)
        save(inventory / 'receipt.json', dict(status='complete_provisional_reference_comparison_inventory', unique_distinct_model_pairs=5, all_reference_comparison_models=6, artifacts={p.name: sha(p) for p in inventory.iterdir()}))
        ir = root / 'inventory-proof.json'; save(ir, dict(status='passed_full_reference_comparison_ledger_and_native_model_readback', producer_receipt_sha256=sha(inventory / 'receipt.json'), distinct_model_pairs=5, unique_models=6))
        input_sources = {}
        for label, subset in [('primary', ['A', 'B', 'C']), ('additional', ['D', 'E', 'F'])]:
            folder = root / label; folder.mkdir(); records = []
            for name in subset:
                for mask in ['full', 'plddt70']:
                    length = lengths[name]; retained = length if mask == 'full' else min(80, length)
                    if name in ['A', 'B'] and mask == 'plddt70': retained = 50
                    if name == 'D' and mask == 'plddt70': retained = 2
                    records.append(dict(model_id=name, version=versions[name], mask=mask, original_length=length,
                                        retained_residues=retained, sequence='A' * retained, original_positions=list(range(1, retained + 1)),
                                        source_sha256=hashlib.sha256(name.encode()).hexdigest(), status='ready' if retained >= 3 else 'too_few_retained_residues'))
            (folder / 'inputs.jsonl').write_text(''.join(json.dumps(r) + '\n' for r in records))
            inp = root / (label + '-plan.json'); save(inp, dict(output=str(folder)))
            save(folder / 'receipt.json', dict(status='complete_duplication_alignment_input_materialization' if label == 'primary' else 'complete_additional_reference_alignment_input_materialization',
                 plan_sha256=sha(inp), models=3, input_dispositions=6, counts=dict(Counter(r['mask'] + ':' + r['status'] for r in records)), artifacts={'inputs.jsonl': sha(folder / 'inputs.jsonl')}))
            proof = root / (label + '-proof.json'); save(proof, dict(status='passed_full_duplication_alignment_input_readback' if label == 'primary' else 'passed_full_reference_alignment_input_readback',
                 **{('source_receipt_sha256' if label == 'primary' else 'producer_receipt_sha256'): sha(folder / 'receipt.json')}))
            input_sources[label] = dict(inputs=str(folder), input_plan=str(inp), readback=str(proof))
        native_plan = root / 'native-plan.json'; save(native_plan, dict(inventory=str(inventory), inventory_readback=str(ir), input_sources=input_sources, pins={}))
        union_plan = root / 'union-plan.json'; save(union_plan, dict(native_plan=str(native_plan), output=str(union)))
        states = []
        for ix, pair in enumerate(pairs):
            ends = [(pair['model_a'], pair['version_a']), (pair['model_b'], pair['version_b'])]
            for mask in ['full', 'plddt70']:
                for order in [0, 1]:
                    directed = ends if order == 0 else ends[::-1]; source_order = 1 - order if ix == 0 else order
                    status = 'aligned'; count = [100, 71, 101, 60, 91][ix] if mask == 'full' else [50, 50, 50, 2, 70][ix]
                    if mask == 'plddt70':
                        if ix == 1 and order == 1: status = 'timeout'
                        if ix == 2: status = 'native_error' if order == 0 else 'parse_error'
                        if ix == 3: status = 'input_unavailable'
                    why = [] if status == 'aligned' else [status]
                    n = g = metrics = None
                    if status == 'aligned':
                        n = dict(pair_key=pair['pair_key'], mask=mask, order=str(source_order), aligned_length=str(count), rmsd_status='within_printed_rounding',
                                 rmsd_recomputed=str(.1 + .01 * order), sequence_identity_exact='0.5', joint_plddt70_fraction='0.9',
                                 tm_left_native='0.6', tm_right_native='0.8', coverage_left=str(count / (lengths[directed[0][0]] if mask == 'full' else 50)),
                                 coverage_right=str(count / (lengths[directed[1][0]] if mask == 'full' else 80)))
                        # Masked retained lengths differ by model, including current endpoint reversal.
                        for side, end in zip(['left', 'right'], directed):
                            retained = lengths[end[0]] if mask == 'full' else 50 if end[0] in ['A', 'B'] else 80
                            n['coverage_' + side] = str(count / retained)
                        g = dict(pair_key=pair['pair_key'], mask=mask, order=str(source_order), aligned_length=str(count), rmsd_status=n['rmsd_status'], geometry_status='unique_at_numeric_tolerance')
                        metrics = dict(aligned_length=count)
                        if ix == 1 and mask == 'full' and order == 0:
                            n['aligned_length'] = g['aligned_length'] = '2'; metrics['aligned_length'] = 2
                            n['rmsd_status'] = g['rmsd_status'] = 'outside_printed_rounding'; g['geometry_status'] = 'degenerate_at_numeric_tolerance'
                            why = ['rmsd_discrepancy', 'fewer_than_three_pairs', 'nonunique_rotation']
                    states.append(dict(pair_key=pair['pair_key'], model_a=ends[0][0], version_a=ends[0][1], model_b=ends[1][0], version_b=ends[1][1], mask=mask, order=order,
                         directed_endpoints=[dict(model_id=m, version=v) for m, v in directed], source_order=source_order, source_native_status=status,
                         selected_source='primary_expanded_completed' if ix == 0 else 'reference_old' if ix == 1 else 'reference_new_native',
                         source_checkpoint=f'synthetic/{pair["pair_key"]}-{mask}-{source_order}.json', source_checkpoint_sha256='a' * 64,
                         source_numeric=n, source_geometry=g, native_metrics=metrics, numerical_exclusion_reasons=why, numerical_usable=not why))
        states.sort(key=lambda r: (r['pair_key'], r['mask'], r['order']))
        union_table = union / 'reference_measurement_dispositions.jsonl'
        union_table.write_text(''.join(json.dumps(r) + '\n' for r in states))
        summary = dict(full_pairs=5, directed_dispositions=20, source_dispositions=dict(Counter(r['selected_source'] for r in states)),
                       native_status_counts=dict(Counter(r['mask'] + ':' + r['source_native_status'] for r in states)),
                       numerical_counts=dict(Counter('usable' if r['numerical_usable'] else ';'.join(r['numerical_exclusion_reasons']) for r in states)))
        save(union / 'receipt.json', dict(status='complete_full_reference_measurement_union_pending_independent_readback', plan_sha256=sha(union_plan), **summary,
                                        artifacts={union_table.name: sha(union_table)}))
        up = root / 'union-proof.json'; save(up, dict(status='passed_full_reference_measurement_union_readback', plan_sha256=sha(union_plan), producer_receipt_sha256=sha(union / 'receipt.json'), **summary))
        archive = root / 'union-archive.json'; save(archive, dict(services=[{'synthetic_stub': True}] * 2, summary=summary, source_hashes={str(union_table): sha(union_table)}))
        closed = root / 'union-closed.json'; save(closed, dict(status='complete_verified_full_reference_measurement_union', **summary, full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive),
             producer_receipt=str(union / 'receipt.json'), producer_receipt_sha256=sha(union / 'receipt.json'), independent_readback=str(up), independent_readback_sha256=sha(up), exact_process_journals_checked=2))
        screens = json.loads(Path('metadata/whole_protein_common_fits_plan_20260927.json').read_text())['screens']
        sp = root / 'screens.json'; save(sp, dict(screens=screens))
        plan = root / 'plan.json'; out = root / 'output'
        config = dict(union_plan=str(union_plan), union_completion=str(closed), screens_plan=str(sp), screens=screens, output=str(out), pins={}, expected=dict(pairs=5, models=6, all_input_dispositions=12))
        save(plan, config)
        run([sys.executable, 'scripts/screen_full_reference_orders_and_coverage.py', '--plan', str(plan)])
        read_command = [sys.executable, 'scripts/readback_full_reference_orders_and_coverage.py', '--plan', str(plan), '--output']
        run(read_command + [str(root / 'proof.json')])
        original = list(csv.DictReader((out / 'pair_mask_order_coverage.tsv').open(), delimiter='\t'))
        lookup = {(r['pair_key'], r['mask']): r for r in original}
        ab = lookup[pairs[0]['pair_key'], 'full']; assert ab['n50_c50_pass'] == '1' and ab['order0_source_order'] == '1' and ab['order0_tm_a'] == '0.6'
        ab70 = lookup[pairs[0]['pair_key'], 'plddt70']; assert ab70['order0_retained_input_coverage_a'] == '1.0' and ab70['order0_original_coverage_a'] == '0.25' and ab70['n30_c50_pass'] == '0'
        ac = lookup[pairs[1]['pair_key'], 'full']; assert ac['order0_native_aligned_length'] == '2' and ac['order0_aligned_length'] == '' and ac['n30_c50_exclusions'] == 'order0_not_numerically_usable'
        ef = lookup[pairs[4]['pair_key'], 'plddt70']; assert ef['n50_c70_pass'] == '0' and ef['n50_c50_pass'] == '1'
        assert lookup[pairs[4]['pair_key'], 'full']['n50_c90_pass'] == '1'
        original_receipt = json.loads((out / 'receipt.json').read_text()); rejected = []
        for label in ['retained_denominator_as_original', 'source_order_as_current_direction', 'favorable_order_pass', 'cleared_quarantine', 'changed_version', 'deleted_unavailable_row', 'changed_row_hash']:
            rows = [dict(r) for r in original]
            if label == 'retained_denominator_as_original':
                row = next(r for r in rows if r['pair_key'] == pairs[0]['pair_key'] and r['mask'] == 'plddt70'); row['order0_original_coverage_a'] = row['order0_retained_input_coverage_a']
            elif label == 'source_order_as_current_direction':
                row = next(r for r in rows if r['pair_key'] == pairs[0]['pair_key']); row['order0_tm_a'], row['order0_tm_b'] = row['order0_tm_b'], row['order0_tm_a']
            elif label in ['favorable_order_pass', 'cleared_quarantine']:
                row = next(r for r in rows if r['pair_key'] == pairs[1]['pair_key'] and r['mask'] == 'full')
                if label == 'favorable_order_pass': row['n50_c70_pass'] = '1'; row['n50_c70_exclusions'] = ''
                else: row['order0_numerical_exclusion_reasons'] = ''; row['order0_status'] = 'aligned'
            elif label == 'changed_version': rows[0]['version_a'] = '99'
            elif label == 'deleted_unavailable_row': rows = [r for r in rows if not (r['pair_key'] == pairs[3]['pair_key'] and r['mask'] == 'plddt70')]
            else: rows[0]['order0_union_row_sha256'] = '0' * 64
            table(out / 'pair_mask_order_coverage.tsv', rows)
            receipt = dict(original_receipt); receipt['artifacts'] = {'pair_mask_order_coverage.tsv': sha(out / 'pair_mask_order_coverage.tsv')}; save(out / 'receipt.json', receipt)
            result = subprocess.run(read_command + [str(root / (label + '.json'))], capture_output=True, text=True)
            assert result.returncode != 0 and not (root / (label + '.json')).exists(), label; rejected.append(label)
        # Rehashed incorrect source original-length manifests must also fail before completion.
        source_file = Path(input_sources['primary']['inputs']) / 'inputs.jsonl'
        source_rows = [json.loads(line) for line in source_file.read_text().splitlines()]; source_rows[0]['original_length'] = 99
        source_file.write_text(''.join(json.dumps(r) + '\n' for r in source_rows))
        rp = source_file.parent / 'receipt.json'; r = json.loads(rp.read_text()); r['artifacts']['inputs.jsonl'] = sha(source_file); save(rp, r)
        proof = Path(input_sources['primary']['readback']); pp = json.loads(proof.read_text()); pp['source_receipt_sha256'] = sha(rp); save(proof, pp)
        bad_plan = root / 'bad-source-plan.json'; save(bad_plan, dict(config, output=str(root / 'bad-source-output')))
        bad = subprocess.run([sys.executable, 'scripts/screen_full_reference_orders_and_coverage.py', '--plan', str(bad_plan)], capture_output=True, text=True)
        assert bad.returncode != 0 and not (root / 'bad-source-output' / 'receipt.json').exists()
        print(json.dumps(dict(status='passed_full_synthetic_reference_order_and_original_coverage_checks', pairs=5, directed_states=20, pair_mask_rows=10, screen_decisions=60,
             rejected_rehashed_exports=rejected, rejected_rehashed_source='incorrect_original_length',
             scope='Synthetic metric/source/proof/journal fixture data only, not production qualification or a biological pilot. Both masks/orders, unequal versions, nonlexical current endpoints, reversed old source directions, exact integer/decimal coverage boundaries, retained-denominator trap, all five native statuses, all three numerical exclusions and blank excluded metrics exercised. Producer and independent decimal-ceiling reader both passed.'), indent=2))


if __name__ == '__main__': main()
