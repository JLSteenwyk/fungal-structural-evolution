#!/usr/bin/env python3
"""Check literal geometry, score boundaries, exact serialization and corruption dispositions."""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile

import numpy as np

from ancestral_chain_attempt import sha
from audit_full_atlas_coordinates_v1 import profile, shard


CIF = """data_literal
_entity_poly.pdbx_seq_one_letter_code_can 'AG'
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_seq_id
_atom_site.label_asym_id
_atom_site.label_alt_id
_atom_site.pdbx_PDB_model_num
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.type_symbol
ATOM N ALA 1 A . 1 0.0001 2.0 3.0 1.0 45.0 N
ATOM CA ALA 1 A . 1 1.0001 2.0 3.0 1.0 49.99 C
ATOM C ALA 1 A . 1 2.0001 2.0 3.0 1.0 45.0 C
ATOM O ALA 1 A . 1 3.0001 2.0 3.0 1.0 45.0 O
ATOM CA GLY 2 A . 1 4.0001 5.0 6.0 1.0 70.0 C
#
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    checks = []
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        path = root/'literal.cif'
        path.write_text(CIF)
        model = dict(path=str(path), sha256=sha(path), sequence_sha256=hashlib.sha256(b'AG').hexdigest(),
                     length=2, mean_ca_plddt=59.995, model_id='literal', version=1)
        blob, summary = profile(model)
        with np.load(io.BytesIO(blob), allow_pickle=False) as arrays:
            assert arrays['ca_xyz'].dtype == np.dtype('float64')
            assert np.array_equal(arrays['ca_xyz'], [[1.0001, 2., 3.], [4.0001, 5., 6.]])
            assert np.array_equal(arrays['ca_plddt'], [49.99, 70.])
            assert arrays['sequence_ascii'].tobytes() == b'AG'
            assert np.array_equal(arrays['missing_backbone'], [False, True])
        assert summary['plddt_ge50'] == 1 and summary['plddt_ge70'] == 1 and summary['plddt_ge90'] == 0
        assert summary['missing_backbone_residues'] == 1
        assert abs(summary['mean_coordinate_ca_plddt']-59.995) < 1e-12
        checks += ['literal_C_alpha_geometry_exact_float64_without_PDB_rounding',
                   'literal_sequence_and_inclusive50_70_90_confidence_boundaries',
                   'missing_backbone_retained_without_silent_drop']
        changes = [('nonfinite_coordinate', CIF.replace('4.0001 5.0', 'nan 5.0')),
                   ('out_of_bounds_confidence', CIF.replace('70.0 C', '100.01 C')),
                   ('residue_sequence_mismatch', CIF.replace('CA GLY 2', 'CA ALA 2')),
                   ('incomplete_C_alpha', CIF.replace('ATOM CA GLY', 'ATOM N GLY')),
                   ('multiple_chains', CIF.replace('GLY 2 A', 'GLY 2 B'))]
        for label, text in changes:
            changed = root/(label+'.cif'); changed.write_text(text)
            bad = dict(model, path=str(changed), sha256=sha(changed))
            try:
                profile(bad)
            except ValueError:
                checks.append(label+'_rejected')
            else:
                raise AssertionError(label+' accepted')
        job = root/'job.jsonl'
        job.write_text(json.dumps(dict(source='ESMFold', model=model))+'\n')
        out = root/'shards'; out.mkdir()
        result = shard(job, out, 0)
        assert result['counts']['ESMFold']['valid_models'] == 1
        assert result['counts']['ESMFold']['valid_residues'] == 2
        rows = [json.loads(line) for line in (out/'job.jsonl').read_text().splitlines()]
        assert len(rows) == 1 and rows[0]['source_sha256'] == model['sha256']
        assert rows[0]['pae_evaluated'] is False
        with tarfile.open(out/'job.tar') as archive:
            assert archive.getnames() == ['ESMFold-'+model['sequence_sha256']+'.npz']
            assert archive.extractfile(rows[0]['member']).read() == blob
        checks.append('source_bound_tar_member_roundtrip_and_full_shard_dispositions')
        corrupt = root/'corrupt.jsonl'
        corrupt.write_text(json.dumps(dict(source='AFDB', model=dict(model, sha256='0'*64)))+'\n')
        try:
            shard(corrupt, out, 0)
        except ValueError as error:
            assert str(error) == 'Source coordinate checksum mismatch'
            assert not (out/'corrupt.receipt.json').exists()
            checks.append('immutable_source_corruption_fails_stage_without_completed_shard')
        else:
            raise AssertionError('Changed source bytes accepted')
    pins = {str(path): sha(path) for path in (Path(__file__), Path('scripts/audit_full_atlas_coordinates_v1.py'),
                                             Path('scripts/extract_domain_coordinates.py'))}
    result = dict(status='passed_literal_full_atlas_coordinate_confidence_profile_checks',
                  checked_utc=datetime.now(timezone.utc).isoformat(), checks=checks, source_hashes=pins,
                  scientific_eligibility=False, gpu=False, new_predictions=0,
                  scope='Literal coordinate/confidence/sequence expectations and controlled rejection/corruption '
                        'tests; no corpus subset pilot, source-accuracy calibration or biological qualification. '
                        'Full source-model native audit and independent readback remain required.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
