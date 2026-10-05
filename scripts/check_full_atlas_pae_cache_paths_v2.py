#!/usr/bin/env python3
"""Offline literal regression controls for the unchanged PAE retriever contract."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import verify
from retrieve_marker_pae import ROOT, retrieve, validate_pae


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    pins = {str(p): sha(p) for p in (
        Path(__file__), Path('scripts/retrieve_full_atlas_missing_pae_v1.py'),
        Path('scripts/retrieve_full_atlas_missing_pae_v2.py'),
        Path('scripts/retrieve_marker_pae.py'), Path('scripts/retrieve_matched_models.py'))}
    checks = []
    raw = json.dumps([{'predicted_aligned_error': [[0.25, 4], [12.5, 0.25]],
                       'max_predicted_aligned_error': 12.5}]).encode()
    model = dict(model_id='AF-OFFLINE-LITERAL-F1', version=6, length=2,
                 sequence_sha256=hashlib.sha256(b'AC').hexdigest(),
                 pae_url='https://alphafold.ebi.ac.uk/files/AF-OFFLINE-LITERAL-F1-predicted_aligned_error_v6.json')

    class Response:
        headers = {}
        def __enter__(self): return self
        def __exit__(self, *args): return False
        def read(self): return raw

    with tempfile.TemporaryDirectory(prefix='offline-pae-paths-', dir=ROOT/'results') as directory:
        absolute = Path(directory)
        relative = absolute.relative_to(ROOT)
        with patch('retrieve_marker_pae.urllib.request.urlopen', return_value=Response()) as http:
            try:
                retrieve(model, relative)
            except ValueError as error:
                assert 'one path is relative and the other is absolute' in str(error)
                checks.append('original_relative_path_receipt_failure_reproduced')
            else:
                raise AssertionError('Original relative path did not fail')
            assert http.call_count == 1
            assert len(list(absolute.glob('*.json.gz'))) == 1
            assert not list(absolute.glob('*.receipt.json'))
            checks.append('original_failure_occurs_after_matrix_write_before_receipt')
            record = retrieve(model, relative.resolve())
            assert record['status'] == 'verified' and http.call_count == 2
            assert record['path'] == str(relative/(model['model_id']+'-v6.json.gz'))
            assert record['json_sha256'] == hashlib.sha256(raw).hexdigest()
            assert all(record[k] == model[k] for k in ('model_id', 'version', 'length', 'sequence_sha256'))
            receipt = absolute/(model['model_id']+'-v6.receipt.json')
            assert json.loads(receipt.read_text()) == record
            checks.append('resolved_cache_preserves_exact_original_relative_receipt_and_identity')
            assert retrieve(model, absolute) == record and http.call_count == 2
            checks.append('verified_cache_readback_avoids_additional_http')
            try:
                retrieve(dict(model, sequence_sha256='0'*64), absolute)
            except ValueError as error:
                assert str(error) == 'Cached PAE provenance mismatch'
                checks.append('changed_sequence_provenance_rejected')
            else:
                raise AssertionError('Changed provenance accepted')
            assert http.call_count == 2
        matrix = validate_pae(raw, 2)
        assert matrix.tolist() == [[0.25, 4], [12.5, 0.25]] and matrix.dtype.str == '<f8'
        assert matrix[0, 1] != matrix[1, 0]
        checks.append('directional_values_and_original_float64_preserved')
        for label, payload, length in (
            ('wrong_dimensions_rejected', raw, 3),
            ('negative_values_rejected', json.dumps([{'predicted_aligned_error': [[0, -1], [2, 0]], 'max_predicted_aligned_error': 2}]).encode(), 2),
            ('declared_maximum_bound_unchanged', json.dumps([{'predicted_aligned_error': [[0, 3], [2, 0]], 'max_predicted_aligned_error': 2}]).encode(), 2)):
            try:
                validate_pae(payload, length)
            except ValueError:
                checks.append(label)
            else:
                raise AssertionError(label)
    verify(pins)
    result = dict(status='passed_offline_literal_pae_cache_path_regression', checks=checks,
                  fixture_http_calls=2, live_http_calls=0, source_hashes=pins,
                  scientific_eligibility=False, gpu=False, new_predictions=0,
                  scope='Literal synthetic response only. Reproduces V1 path error and verifies resolved-cache receipt writing, cached readback and unchanged provenance/shape/directional/max-bound checks. No corpus pilot or live download.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
