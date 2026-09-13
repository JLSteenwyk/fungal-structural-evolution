#!/usr/bin/env python3
"""Combine disjoint frozen local model inventories with matching inference settings."""
import argparse
import json
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohort', action='append', required=True,
                        help='Completed mapping directory=prediction configuration file')
    parser.add_argument('--output', type=Path, required=True)
    a = parser.parse_args()
    if a.output.exists():
        raise FileExistsError('Use an immutable combined inventory')
    records = []; sources = []; sequence_ids = set(); model_ids = set()
    shared = None
    for argument in a.cohort:
        mapping_name, config_name = argument.split('=', 1)
        mapping, config_path = Path(mapping_name), Path(config_name)
        receipt = checked_receipt(mapping)
        if receipt['source_policy']['provider'] != 'local' or receipt['source_policy']['tool'] != 'ESMFold v1':
            raise ValueError('Only local ESMFold snapshots can be combined')
        config = json.loads(config_path.read_text())
        settings = {k: v for k, v in config.items() if k not in ['input_receipt_sha256', 'visible_gpu']}
        if shared is None:
            shared = settings
        elif shared != settings:
            raise ValueError('Inference settings differ')
        config_sha = sha(config_path)
        provenance = json.loads((mapping / 'model_provenance.json').read_text())
        expected = {m['model_id']: m for m in provenance}
        if len(expected) != len(provenance) or len(expected) != receipt['distinct_models']:
            raise ValueError('Mapped model identities differ')
        inventory_path = mapping / 'source_inventory.jsonl'
        if sha(inventory_path) != receipt['source_inventory_sha256']:
            raise ValueError('Frozen inventory changed')
        cohort = [json.loads(line) for line in inventory_path.read_text().splitlines()]
        found = set()
        for row in cohort:
            if row['status'] != 'verified' or len(row['models']) != 1:
                raise ValueError('Expected one verified model per sequence')
            m = row['models'][0]; sid = row['record_id']; mid = m['model_id']
            if sid in sequence_ids or mid in model_ids:
                raise ValueError('Overlapping inventories require explicit reconciliation')
            if mid not in expected or sid != 'S' + m['sequence_sha256']:
                raise ValueError('Inventory/mapping identity differs')
            old = expected[mid]
            for field in ['sequence_sha256', 'sha256', 'path', 'prediction_config_sha256',
                          'prediction_receipt_sha256', 'local_pae_npz_sha256']:
                if m[field] != old[field]:
                    raise ValueError('Model provenance differs: ' + field)
            if m['provider'] != 'local' or m['tool'] != 'ESMFold v1' or m['prediction_config_sha256'] != config_sha:
                raise ValueError('Prediction method/configuration differs')
            if sha(ROOT / m['path']) != m['sha256']:
                raise ValueError('Coordinate file changed')
            sequence_ids.add(sid); model_ids.add(mid); found.add(mid); records.append(row)
        if found != set(expected):
            raise ValueError('Inventory includes a different model set')
        sources.append({'mapping': str(mapping), 'mapping_receipt_sha256': sha(mapping / 'receipt.json'),
                        'inventory_sha256': sha(inventory_path), 'models': len(cohort),
                        'prediction_config': str(config_path), 'prediction_config_sha256': config_sha,
                        'visible_gpu': config.get('visible_gpu'),
                        'input_receipt_sha256': config['input_receipt_sha256']})
    a.output.mkdir(parents=True)
    target = a.output / 'inventory.jsonl'
    with target.open('w') as handle:
        for row in sorted(records, key=lambda x: x['record_id']):
            handle.write(json.dumps(row, sort_keys=True) + '\n')
    result = {'status': 'complete_disjoint_same_method_inventory_union', 'models': len(records),
              'source_cohorts': sources, 'shared_inference_settings': shared,
              'script_sha256': sha(Path(__file__)), 'artifacts': {'inventory.jsonl': sha(target)},
              'scope': 'All frozen model identities, provenance and coordinate hashes verified; no overlapping sequences or models. Only input receipt and physical GPU identity differ between configurations. No coordinate averaging or prediction rerun. Remapping, combined qualified encodings and paired inference inputs remain separate; matching settings do not prove absence of batch effects.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
