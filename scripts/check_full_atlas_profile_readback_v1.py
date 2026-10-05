#!/usr/bin/env python3
"""Exercise independent profile readback against literal coordinates and exact 50/90 cutoffs."""
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
from check_full_atlas_coordinate_profiles_v1 import CIF
from readback_full_atlas_coordinate_profiles_v1 import independent_arrays, check_profile, read_shard


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); assert not args.receipt.exists()
    checks = []
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        original = root/'literal.cif'
        original.write_text(CIF.replace('49.99 C', '50.0 C').replace('70.0 C', '90.0 C'))
        model = dict(path=str(original), sha256=sha(original), sequence_sha256=hashlib.sha256(b'AG').hexdigest(),
                     length=2, mean_ca_plddt=70., model_id='literal', version=1)
        expected = dict(sequence_ascii=np.asarray([65, 71], dtype=np.uint8),
                        ca_xyz=np.asarray([[1.0001, 2., 3.], [4.0001, 5., 6.]], dtype=np.float64),
                        ca_plddt=np.asarray([50., 90.], dtype=np.float64),
                        missing_backbone=np.asarray([False, True], dtype=np.bool_))
        decoded = independent_arrays(model)
        assert set(decoded) == set(expected) and all(np.array_equal(decoded[k], v) for k, v in expected.items())
        buffer = io.BytesIO(); np.savez_compressed(buffer, **expected); blob = buffer.getvalue()
        row = dict(source='ESMFold', sequence_sha256=model['sequence_sha256'], model_id='literal', version=1,
                   source_path=model['path'], source_sha256=model['sha256'], length=2, prediction_config_sha256='',
                   confidence_representation='Original mmCIF C-alpha B_iso_or_equiv; score, not calibrated accuracy',
                   pae_evaluated=False, status='coordinate_and_residue_profile_exported',
                   member='ESMFold-'+model['sequence_sha256']+'.npz', npz_sha256=hashlib.sha256(blob).hexdigest(),
                   residues=2, missing_backbone_residues=1, mean_coordinate_ca_plddt=70., minimum_ca_plddt=50.,
                   maximum_ca_plddt=90., catalog_mean_ca_plddt=70., coordinate_minus_catalog_mean_ca_plddt=0.,
                   plddt_ge50=2, plddt_ge70=1, plddt_ge90=1,
                   array_sha256={k: hashlib.sha256(v.dtype.str.encode()+str(v.shape).encode()+v.tobytes()).hexdigest()
                                 for k, v in expected.items()})
        check_profile(blob, row, decoded, model)
        checks += ['independent_literal_source_coordinates_sequence_and_backbone',
                   'exact50_and90_inclusive_counts_with_independent_scalar_recount']
        wrong = dict(row, plddt_ge90=0)
        try:
            check_profile(blob, wrong, decoded, model)
        except AssertionError:
            checks.append('wrong_threshold_count_rejected')
        else:
            raise AssertionError('Incorrect 90 cutoff count accepted')
        altered = {k: v.copy() for k, v in expected.items()}; altered['ca_xyz'][0, 0] += .0001
        buffer = io.BytesIO(); np.savez_compressed(buffer, **altered); changed = buffer.getvalue()
        changed_row = dict(row, npz_sha256=hashlib.sha256(changed).hexdigest(),
                           array_sha256={k: hashlib.sha256(v.dtype.str.encode()+str(v.shape).encode()+v.tobytes()).hexdigest()
                                         for k, v in altered.items()})
        try:
            check_profile(changed, changed_row, decoded, model)
        except AssertionError:
            checks.append('self_consistent_modified_array_rejected_against_original_source')
        else:
            raise AssertionError('Changed source geometry accepted')
        try:
            independent_arrays(dict(model, sha256='0'*64))
        except RuntimeError:
            checks.append('immutable_source_corruption_is_hard_failure')
        else:
            raise AssertionError('Corrupt source accepted')
        job = root/'shard_00000.jsonl'; job.write_text(json.dumps(dict(source='ESMFold', model=model))+'\n')
        shards = root/'shards'; shards.mkdir()
        with tarfile.open(shards/'shard_00000.tar', 'w', format=tarfile.USTAR_FORMAT) as archive:
            member = tarfile.TarInfo(row['member']); member.size=len(blob); member.mtime=0
            archive.addfile(member, io.BytesIO(blob))
        rows = shards/'shard_00000.jsonl'; rows.write_text(json.dumps(row)+'\n')
        zero = dict(models=0, valid_models=0, rejected_models=0, source_coordinate_bytes=0,
                    valid_residues=0, missing_backbone_residues=0, plddt_ge50=0, plddt_ge70=0, plddt_ge90=0)
        counts = dict(AFDB=zero, ESMFold=dict(models=1, valid_models=1, rejected_models=0,
                      source_coordinate_bytes=original.stat().st_size, valid_residues=2,
                      missing_backbone_residues=1, plddt_ge50=2, plddt_ge70=1, plddt_ge90=1))
        proof = dict(job=str(job), receipt=dict(job_sha256=sha(job), counts=counts,
                     artifacts={p.name: sha(p) for p in (rows, shards/'shard_00000.tar')}))
        assert read_shard(proof, shards) == counts
        checks.append('complete_literal_shard_identity_membership_and_counts')
    paths = [Path(__file__), Path('scripts/readback_full_atlas_coordinate_profiles_v1.py'),
             Path('scripts/check_full_atlas_coordinate_profiles_v1.py')]
    result = dict(status='passed_independent_literal_full_atlas_profile_readback_checks',
                  checked_utc=datetime.now(timezone.utc).isoformat(), checks=checks,
                  source_hashes={str(p): sha(p) for p in paths}, scientific_eligibility=False,
                  gpu=False, new_predictions=0,
                  scope='Literal-source reader controls, exact50/90 boundaries, full tar/job/profile '
                        'readback and self-consistent geometry/count/source-corruption rejection. '
                        'No corpus pilot or full native audit/readback completion. Shared CIF lexical '
                        'parser; independent source reconstruction and scalar counts.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
