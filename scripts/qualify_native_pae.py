#!/usr/bin/env python3
"""Add validated PAE to a completed coordinate audit without repeating CIF parsing."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import numpy as np
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT, sha
from retrieve_marker_pae import validate_pae


def context_maximum(valid, partner_1based, pae):
    i = np.flatnonzero(valid)
    j = partner_1based[i] - 1
    if np.any(i < 1) or np.any(i >= len(valid)-1) or np.any(j < 1) or np.any(j >= len(valid)-1):
        raise ValueError('Invalid six-residue context')
    context = np.stack([i-1, i, i+1, j-1, j, j+1], axis=1)
    maximum = np.full(len(valid), np.nan)
    maximum[i] = pae[context[:, :, None], context[:, None, :]].max(axis=(1, 2))
    return maximum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coordinates', type=Path, required=True)
    parser.add_argument('--pae', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable PAE-qualified output')
    coordinate = checked_receipt(args.coordinates)
    confidence = checked_receipt(args.pae)
    if coordinate['status'] != 'complete_native_3di_coordinate_audit' or coordinate['pae_receipt_sha256'] is not None:
        raise ValueError('Completed coordinate-only audit required')
    if coordinate['mapping_receipt_sha256'] != confidence['mapping_receipt_sha256']:
        raise ValueError('Coordinate and PAE mapping snapshots differ')
    native = Path(coordinate['native_path'])
    if sha(native / 'config.json') != coordinate['native_config_sha256']:
        raise ValueError('Changed native configuration')
    for name, checksum in coordinate['native_artifacts'].items():
        if sha(native / name) != checksum:
            raise ValueError('Changed native artifact')
    with (args.coordinates / 'model_summary.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    models = json.loads((native / 'model_provenance.json').read_text())
    expected = {(r['model_id'], str(r['version'])): r for r in models}
    pae_rows = json.loads((args.pae / 'pae_manifest.json').read_text())
    paes = {(r['model_id'], str(r['version'])): r for r in pae_rows}
    keys = [(r['model_id'], str(r['version'])) for r in rows]
    if len(expected) != len(models) or len(paes) != len(pae_rows) or len(set(keys)) != len(rows) or set(keys) != set(expected) or set(paes) != set(expected):
        raise ValueError('Model universe differs or repeats an identity')
    args.output = args.output.resolve()
    args.output.mkdir(parents=True)
    summaries = []
    integer_fields = ['length', 'valid_states', 'invalid_states', 'valid_focal_plddt70', 'valid_feature_plddt70']
    for row in rows:
        key = (row['model_id'], str(row['version']))
        model, pae_row = expected[key], paes[key]
        if row['sequence_sha256'] != model['sequence_sha256'] or pae_row['sequence_sha256'] != model['sequence_sha256'] or pae_row['status'] != 'verified':
            raise ValueError('Model sequence or confidence status differs')
        path = ROOT / row['encoding_path']
        if sha(path) != row['encoding_sha256']:
            raise ValueError('Changed audited coordinate encoding')
        with np.load(path, allow_pickle=False) as data:
            arrays = {name: data[name].copy() for name in data.files}
        if 'feature_max_pae' in arrays:
            raise ValueError('Unexpected pre-existing PAE qualification')
        length = int(row['length'])
        pae_path = ROOT / pae_row['path']
        if sha(pae_path) != pae_row['gzip_sha256']:
            raise ValueError('Changed PAE artifact')
        pae = validate_pae(gzip.decompress(pae_path.read_bytes()), length)
        maximum = context_maximum(arrays['valid'], arrays['partner_residue_1based'], pae)
        arrays['feature_max_pae'] = maximum
        output = args.output / (row['model_name'] + '.npz')
        np.savez_compressed(output, **arrays)
        summary = dict(row)
        summary.update({k: int(row[k]) for k in integer_fields})
        summary.update(valid_feature_plddt70_pae10=int((arrays['valid'] & (arrays['feature_min_plddt'] >= 70) & (maximum <= 10)).sum()),
                       encoding_path=str(output.relative_to(ROOT)), encoding_sha256=sha(output))
        summaries.append(summary)
    table = args.output / 'model_summary.tsv'
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, list(summaries[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(summaries)
    result = dict(coordinate)
    result.update(status='complete_native_3di_feature_audit', confidence_stage='plddt_and_pae',
        totals={k: sum(r[k] for r in summaries) for k in integer_fields + ['valid_feature_plddt70_pae10']},
        coordinate_audit_receipt_sha256=sha(args.coordinates / 'receipt.json'),
        coordinate_audit_path=str(args.coordinates), pae_receipt_sha256=sha(args.pae / 'receipt.json'),
        coordinate_audit_script_sha256=coordinate['script_sha256'], script_sha256=sha(Path(__file__)),
        numpy_version=np.__version__, artifacts={'model_summary.tsv': sha(table)})
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
