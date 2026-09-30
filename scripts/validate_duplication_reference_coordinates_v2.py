#!/usr/bin/env python3
"""Validate the exact reference partition after the stronger native gene/model ledger audit."""
import argparse
import fcntl
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from run_ortholog_pair_guide_comparison import sha
from validate_duplication_coordinates import run_shard


def model_map(path):
    rows = [json.loads(line) for line in Path(path).read_text().splitlines()]
    result = {(row['model_id'], row['version']): row for row in rows}
    if len(result) != len(rows):
        raise ValueError('Duplicate model identity')
    return result


def load_additional_models(plan):
    inventory = Path(plan['inventory'])
    base = Path(plan['base_queue'])
    receipt = json.loads((inventory / 'receipt.json').read_text())
    audit = json.loads(Path(plan['inventory_readback']).read_text())
    if receipt['status'] != 'complete_provisional_reference_comparison_inventory':
        raise ValueError('Incomplete reference inventory')
    if (audit['status'] != 'passed_full_reference_comparison_ledger_and_native_model_readback'
            or audit['producer_receipt_sha256'] != sha(inventory / 'receipt.json')):
        raise ValueError('Reference audit does not bind inventory')
    for name, digest in receipt['artifacts'].items():
        if sha(inventory / name) != digest:
            raise ValueError('Changed reference artifact: ' + name)
    base_receipt = json.loads((base / 'receipt.json').read_text())
    if (base_receipt['status'] != 'complete_reviewed_duplication_model_pair_queue'
            or sha(base / 'models.jsonl') != base_receipt['artifacts']['models.jsonl']):
        raise ValueError('Base queue binding differs')
    for path in [base / 'receipt.json', base / 'models.jsonl']:
        if receipt['source_pins'][str(path)] != sha(path):
            raise ValueError('Reference inventory used another base queue')
    all_models = model_map(inventory / 'models.jsonl')
    original = model_map(base / 'models.jsonl')
    extra = model_map(inventory / 'additional_models.jsonl')
    if len(original) != base_receipt['unique_models']:
        raise ValueError('Base model count differs')
    expected = {key: row for key, row in all_models.items() if key not in original}
    if extra != expected:
        raise ValueError('Additional models are not the exact set difference')
    if any(row != original[key] for key, row in all_models.items() if key in original):
        raise ValueError('Reused model provenance differs')
    if (len(extra) != receipt['additional_models']
            or len(extra) != audit['additional_models']
            or len(all_models) != receipt['all_reference_comparison_models']
            or len(all_models) != audit['unique_models']):
        raise ValueError('Model count differs from reviewed inventory')
    return list(extra.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    plan_hash = sha(args.plan)

    def verify():
        if sha(args.plan) != plan_hash:
            raise ValueError('Changed plan')
        for path, digest in plan['pins'].items():
            if sha(path) != digest:
                raise ValueError('Changed pin: ' + path)

    verify()
    models = load_additional_models(plan)
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'run.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        n = plan['models_per_shard']
        if n < 1 or plan['workers'] < 1:
            raise ValueError('Positive shard size and worker count required')
        jobs = [(i, models[start:start+n], str(out), plan_hash,
                 plan['minimum_free_disk_gib'])
                for i, start in enumerate(range(0, len(models), n))]
        totals = Counter()
        receipts = []
        with ProcessPoolExecutor(max_workers=plan['workers']) as pool:
            for receipt in pool.map(run_shard, jobs):
                receipts.append(receipt)
                totals.update(receipt['counts'])
                state = dict(stage='validating_additional_reference_coordinates',
                             completed_shards=len(receipts), total_shards=len(jobs),
                             counts=dict(totals))
                (out / 'state.json').write_text(json.dumps(state) + '\n')
                print(json.dumps(state), flush=True)
        verify()
        if sum(totals.values()) != len(models):
            raise ValueError('Incomplete model dispositions')
        result = dict(
            status='complete_duplication_coordinate_validation_with_dispositions',
            plan_sha256=plan_hash, models=len(models), counts=dict(totals),
            shards=receipts, script_sha256=sha(__file__),
            inventory_receipt_sha256=sha(Path(plan['inventory']) / 'receipt.json'),
            scope='Only the exact additional-model partition of the provisional '
                  'sister-reference inventory. Same frozen raw-CIF, canonical '
                  'sequence, atom, finite-coordinate and confidence-summary '
                  'checks as the primary duplication validation. Content '
                  'rejections retained; provenance failures stop the job. '
                  'Independent coordinate readback and alignments remain '
                  'separate stages. No reference orthology or asymmetry inference.')
        (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
