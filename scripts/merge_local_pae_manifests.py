#!/usr/bin/env python3
"""Bind unchanged, verified local PAE exports to an audited combined mapping."""
import argparse
import json
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT, sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pae', action='append', type=Path, required=True)
    p.add_argument('--encodings', type=Path, required=True)
    p.add_argument('--snapshot', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable output')
    mapping, encoding = checked_receipt(a.snapshot), checked_receipt(a.encodings)
    if encoding['status'] != 'complete_union_of_previously_qualified_encodings' or encoding['mapping_receipt_sha256'] != sha(a.snapshot / 'receipt.json'):
        raise ValueError('Expected audited qualified-encoding union for this mapping')
    expected = {r['pae_receipt_sha256']: r for r in encoding['source_cohorts']}
    if len(expected) != len(encoding['source_cohorts']):
        raise ValueError('Repeated source PAE receipt')
    provenance = json.loads((a.snapshot / 'model_provenance.json').read_text())
    models = {(m['model_id'], str(m['version'])): m for m in provenance}
    if len(models) != len(provenance):
        raise ValueError('Duplicate mapped model')
    merged = {}; sources = []; seen = set()
    for folder in a.pae:
        receipt = checked_receipt(folder); digest = sha(folder / 'receipt.json')
        if digest not in expected or digest in seen:
            raise ValueError('Unexpected or repeated PAE cohort')
        cohort = expected[digest]
        if receipt['mapping_receipt_sha256'] != cohort['mapping_receipt_sha256'] or receipt['status'] != 'complete_local_mapping_bound_pae_export':
            raise ValueError('Source export/mapping provenance differs')
        rows = json.loads((folder / 'pae_manifest.json').read_text())
        if len(rows) != cohort['models']:
            raise ValueError('Incomplete cohort PAE model grid')
        for row in rows:
            key = row['model_id'], str(row['version'])
            if key in merged or key not in models or row['status'] != 'verified':
                raise ValueError('Duplicate, unexpected or unverified model')
            m = models[key]
            for field in ['sequence_sha256', 'length', 'prediction_config_sha256', 'prediction_receipt_sha256']:
                if row[field] != m[field]:
                    raise ValueError('PAE model identity differs: ' + field)
            if row['source_npz_path'] != m['local_pae_npz_path'] or row['source_npz_sha256'] != m['local_pae_npz_sha256']:
                raise ValueError('PAE prediction source differs')
            path = ROOT / row['path']
            if path.stat().st_size != row['compressed_bytes'] or sha(path) != row['gzip_sha256']:
                raise ValueError('Changed exported matrix')
            if sha(ROOT / row['source_npz_path']) != row['source_npz_sha256']:
                raise ValueError('Changed original prediction array')
            merged[key] = row
        seen.add(digest)
        sources.append({'path': str(folder), 'receipt_sha256': digest, 'models': len(rows),
                        'mapping_receipt_sha256': receipt['mapping_receipt_sha256']})
    if seen != set(expected) or set(merged) != set(models) or len(merged) != mapping['distinct_models']:
        raise ValueError('Incomplete source/model union')
    a.output.mkdir(parents=True)
    target = a.output / 'pae_manifest.json'
    target.write_text(json.dumps([merged[k] for k in sorted(merged)], indent=2) + '\n')
    result = {'status': 'complete_union_of_verified_local_pae_exports', 'models': len(merged),
              'mapping_receipt_sha256': sha(a.snapshot / 'receipt.json'),
              'encoding_union_receipt_sha256': sha(a.encodings / 'receipt.json'),
              'source_cohorts': sources, 'script_sha256': sha(Path(__file__)),
              'artifacts': {'pae_manifest.json': sha(target)},
              'interpretation': 'Unchanged source export rows, matrix hashes and original NPZ hashes checked against combined model provenance and qualified-cohort receipts. Original per-cohort lossless-export and full-context audits remain authoritative. No new PAE computation or confidence calibration.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
