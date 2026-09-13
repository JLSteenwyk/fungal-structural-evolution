#!/usr/bin/env python3
"""Check every qualified feature context directly against original local PAE arrays."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT, sha


def table(path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    result = {(r['model_id'], str(r['version'])): r for r in rows}
    if len(result) != len(rows):
        raise ValueError('Repeated model identity')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['snapshot', 'coordinates', 'qualified', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable readback output')
    receipts = {name: checked_receipt(getattr(args, name))
                for name in ['snapshot', 'coordinates', 'qualified']}
    coordinate, qualified = receipts['coordinates'], receipts['qualified']
    if (qualified['status'] != 'complete_native_3di_feature_audit'
            or coordinate['status'] != 'complete_native_3di_coordinate_audit'
            or qualified['coordinate_audit_receipt_sha256'] != sha(args.coordinates / 'receipt.json')
            or qualified['mapping_receipt_sha256'] != sha(args.snapshot / 'receipt.json')
            or coordinate['mapping_receipt_sha256'] != qualified['mapping_receipt_sha256']):
        raise ValueError('Source receipts do not match')
    models = json.loads((args.snapshot / 'model_provenance.json').read_text())
    model_map = {(m['model_id'], str(m['version'])): m for m in models}
    before = table(args.coordinates / 'model_summary.tsv')
    after = table(args.qualified / 'model_summary.tsv')
    if (len(model_map) != len(models) or set(before) != set(after)
            or set(after) != set(model_map) or len(after) != qualified['models']):
        raise ValueError('Model universes differ')
    counts = dict.fromkeys(['length', 'valid_states', 'invalid_states',
                           'valid_focal_plddt70', 'valid_feature_plddt70',
                           'valid_feature_plddt70_pae10'], 0)
    checked = 0
    for key, row in after.items():
        prior, model = before[key], model_map[key]
        for field, value in prior.items():
            if field not in ['encoding_path', 'encoding_sha256'] and row[field] != value:
                raise ValueError('Coordinate summary changed: ' + field)
        paths = [ROOT / prior['encoding_path'], ROOT / row['encoding_path'],
                 ROOT / model['local_pae_npz_path']]
        hashes = [prior['encoding_sha256'], row['encoding_sha256'], model['local_pae_npz_sha256']]
        if any(sha(p) != h for p, h in zip(paths, hashes)):
            raise ValueError('Changed encoding or original confidence array')
        with np.load(paths[0], allow_pickle=False) as old, np.load(paths[1], allow_pickle=False) as new, np.load(paths[2], allow_pickle=False) as raw:
            if set(new.files) != set(old.files) | {'feature_max_pae'}:
                raise ValueError('Unexpected qualified arrays')
            for name in old.files:
                options = {'equal_nan': True} if old[name].dtype.kind == 'f' else {}
                if not np.array_equal(old[name], new[name], **options):
                    raise ValueError('Coordinate array changed: ' + name)
            if str(raw['sequence']) != str(new['sequence']):
                raise ValueError('Original prediction sequence differs')
            valid = new['valid']; n = len(valid)
            if n != model['length'] or n != int(row['length']):
                raise ValueError('Length differs')
            focal = np.flatnonzero(valid)
            partner = new['partner_residue_1based'][focal] - 1
            if np.any(focal < 1) or np.any(focal >= n-1) or np.any(partner < 1) or np.any(partner >= n-1):
                raise ValueError('Invalid feature neighborhood')
            pae = raw['pae']
            if pae.shape != (n, n) or not np.isfinite(pae).all() or (pae < 0).any():
                raise ValueError('Invalid original PAE')
            # Independent accumulation of all 36 ordered residue pairs; no
            # qualifier helper, exported JSON parser or stacked context tensor.
            maximum = np.zeros(len(focal), dtype=float)
            positions = [focal-1, focal, focal+1, partner-1, partner, partner+1]
            for left in positions:
                for right in positions:
                    maximum = np.maximum(maximum, pae[left, right])
            if (not np.array_equal(maximum, new['feature_max_pae'][focal])
                    or not np.isnan(new['feature_max_pae'][~valid]).all()):
                raise ValueError('PAE context maximum or invalid sentinel differs')
            observed = {'length': n, 'valid_states': int(valid.sum()),
                        'invalid_states': int((~valid).sum()),
                        'valid_focal_plddt70': int((valid & (new['ca_plddt'] >= 70)).sum()),
                        'valid_feature_plddt70': int((valid & (new['feature_min_plddt'] >= 70)).sum()),
                        'valid_feature_plddt70_pae10': int(((new['feature_min_plddt'][focal] >= 70) & (maximum <= 10)).sum())}
            for name, value in observed.items():
                if value != int(row[name]):
                    raise ValueError('Summary count differs: ' + name)
                counts[name] += value
        checked += 1
        if checked % 500 == 0:
            print('Checked', checked, 'of', len(models), flush=True)
    if counts != qualified['totals']:
        raise ValueError('Aggregate totals differ')
    result = {'status': 'passed_full_original_pae_context_readback', 'models': checked,
              'totals': counts, 'source_receipts': {name: sha(getattr(args, name) / 'receipt.json')
              for name in receipts}, 'script_sha256': sha(Path(__file__)),
              'interpretation': 'Every valid six-residue context checked against original local NPZ directional PAE by 36 ordered-pair maxima; all pre-existing coordinate arrays and summary counts preserved. Does not repeat native partner selection or establish confidence calibration or biological accuracy.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
