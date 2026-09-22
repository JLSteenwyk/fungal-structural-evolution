#!/usr/bin/env python3
"""Check every qualified AlphaFold feature context against version-bound PAE JSON."""
import argparse
import csv
import gzip
import hashlib
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
    for name in ['snapshot', 'coordinates', 'qualified', 'pae', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable readback output')
    receipts = {name: checked_receipt(getattr(args, name))
                for name in ['snapshot', 'coordinates', 'qualified', 'pae']}
    coordinate, qualified = receipts['coordinates'], receipts['qualified']
    confidence = receipts['pae']
    if (qualified['status'] != 'complete_native_3di_feature_audit'
            or coordinate['status'] != 'complete_native_3di_coordinate_audit'
            or qualified['coordinate_audit_receipt_sha256'] != sha(args.coordinates / 'receipt.json')
            or qualified['mapping_receipt_sha256'] != sha(args.snapshot / 'receipt.json')
            or coordinate['mapping_receipt_sha256'] != qualified['mapping_receipt_sha256']
            or qualified['pae_receipt_sha256'] != sha(args.pae / 'receipt.json')
            or confidence['mapping_receipt_sha256'] != qualified['mapping_receipt_sha256']
            or confidence['models_failed'] != 0):
        raise ValueError('Source receipts do not match')
    models = json.loads((args.snapshot / 'model_provenance.json').read_text())
    model_map = {(m['model_id'], str(m['version'])): m for m in models}
    pae_rows = json.loads((args.pae / 'pae_manifest.json').read_text())
    paes = {(r['model_id'], str(r['version'])): r for r in pae_rows}
    if len(paes) != len(pae_rows) or set(paes) != set(model_map):
        raise ValueError('PAE model universe differs')
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
        pae_row = paes[key]
        if (pae_row['status'] != 'verified' or pae_row['sequence_sha256'] != model['sequence_sha256']
                or pae_row['length'] != model['length'] or pae_row['url'] != model['pae_url']):
            raise ValueError('PAE provenance differs')
        paths = [ROOT / prior['encoding_path'], ROOT / row['encoding_path'], ROOT / pae_row['path']]
        hashes = [prior['encoding_sha256'], row['encoding_sha256'], pae_row['gzip_sha256']]
        if any(sha(p) != h for p, h in zip(paths, hashes)):
            raise ValueError('Changed encoding or original confidence array')
        raw = gzip.decompress(paths[2].read_bytes())
        if hashlib.sha256(raw).hexdigest() != pae_row['json_sha256']:
            raise ValueError('Uncompressed PAE checksum differs')
        payload = json.loads(raw)
        if not isinstance(payload, list) or len(payload) != 1:
            raise ValueError('Expected one monomer PAE object')
        pae = np.asarray(payload[0]['predicted_aligned_error'], dtype=float)
        declared = float(payload[0]['max_predicted_aligned_error'])
        if not np.isfinite(declared) or declared < 0 or pae.max() > declared + .51:
            raise ValueError('Invalid declared PAE maximum')
        with np.load(paths[0], allow_pickle=False) as old, np.load(paths[1], allow_pickle=False) as new:
            if set(new.files) != set(old.files) | {'feature_max_pae'}:
                raise ValueError('Unexpected qualified arrays')
            for name in old.files:
                options = {'equal_nan': True} if old[name].dtype.kind == 'f' else {}
                if not np.array_equal(old[name], new[name], **options):
                    raise ValueError('Coordinate array changed: ' + name)
            if hashlib.sha256(str(new['sequence']).encode()).hexdigest() != model['sequence_sha256']:
                raise ValueError('Original prediction sequence differs')
            valid = new['valid']; n = len(valid)
            if n != model['length'] or n != int(row['length']):
                raise ValueError('Length differs')
            focal = np.flatnonzero(valid)
            partner = new['partner_residue_1based'][focal] - 1
            if np.any(focal < 1) or np.any(focal >= n-1) or np.any(partner < 1) or np.any(partner >= n-1):
                raise ValueError('Invalid feature neighborhood')
            if pae.shape != (n, n) or not np.isfinite(pae).all() or (pae < 0).any():
                raise ValueError('Invalid original PAE')
            # Independent accumulation of all 36 ordered residue pairs; no
            # qualifier helper or stacked context tensor.
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
    result = {'status': 'passed_full_afdb_pae_context_readback', 'models': checked,
              'totals': counts, 'source_receipts': {name: sha(getattr(args, name) / 'receipt.json')
              for name in receipts}, 'script_sha256': sha(Path(__file__)),
              'interpretation': 'Every valid six-residue context checked against version-bound AFDB JSON PAE by 36 ordered-pair maxima using a separate JSON/array readback; all pre-existing coordinate arrays and summary counts preserved. Model identity is receipt-bound because PAE JSON does not embed the amino-acid sequence. Does not repeat native partner selection or establish confidence calibration or biological accuracy.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
