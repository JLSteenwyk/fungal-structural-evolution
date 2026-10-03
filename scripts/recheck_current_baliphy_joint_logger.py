#!/usr/bin/env python3
"""Read-only replay of the closed V5 startup and available V3 joint samples."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_joint_sampler_qualification_v3 import inspect, SUCCESS
from run_baliphy_reference_preflight import verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    completion_path = Path('metadata/baliphy_joint_logger_preflight_v5_completed_20261003.json')
    completion = json.loads(completion_path.read_text())
    assert completion['status'] == 'complete_verified_full_joint_logger_startup_v5'
    assert completion['validated_startups'] == 1620 and completion['unsuccessful_startups'] == 0
    archive_path = Path(completion['full_hash_archive'])
    assert sha(archive_path) == completion['full_hash_archive_sha256']
    archive = json.loads(archive_path.read_text())
    verify({'pins': archive['source_hashes']})
    assert len(archive['source_hashes']) == completion['bound_source_hashes']
    plan_path = Path('metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json')
    plan = json.loads(plan_path.read_text())
    verify(plan)
    plan_digest = sha(plan_path)
    jobs = json.loads(Path(plan['jobs']).read_text())
    by_id = {job['chain']['chain_id']: job for job in jobs}
    assert len(by_id) == len(jobs) == 1620
    root = Path(plan['output'])
    # Snapshot immutable completed checkpoints; active roles are never read.
    checkpoint_paths = sorted((root / 'chains').glob('*.json'))
    bindings = {str(completion_path): sha(completion_path), str(archive_path): sha(archive_path),
                str(plan_path): plan_digest, str(Path(__file__)): sha(__file__)}
    rows = []
    for path in checkpoint_paths:
        original_hash = sha(path)
        saved = json.loads(path.read_text())
        cid = saved['chain_id']
        reconstructed = inspect(by_id[cid], Path(saved['native_receipt']), plan_digest,
            plan['mapping'], root / 'frames' / cid, allow_export_creation=False)
        assert saved == reconstructed, cid
        assert sha(path) == original_hash
        bindings[str(path)] = original_hash
        receipt_path = Path(saved['native_receipt'])
        bindings[str(receipt_path)] = sha(receipt_path)
        for name, digest in json.loads(receipt_path.read_text())['artifacts'].items():
            artifact = receipt_path.parent / name
            assert sha(artifact) == digest
            bindings[str(artifact)] = digest
        for frame in saved['joint_frames']:
            bindings[frame['projection_array']] = frame['projection_array_sha256']
        rows.append(saved)
    verify({'pins': bindings})
    status_counts = dict(Counter(row['status'] for row in rows))
    frames = [frame for row in rows for frame in row['joint_frames']]
    result = dict(status='completed_read_only_current_joint_logger_error_recheck',
        checked_utc=datetime.now(timezone.utc).isoformat(),
        closed_startups_checked=completion['validated_startups'],
        startup_closure_bindings_rehashed=len(archive['source_hashes']),
        available_completed_roles_replayed=len(rows), expected_sampler_roles=1620,
        status_counts=status_counts,
        unsuccessful_or_invalid_available_roles=sum(row['status'] != SUCCESS for row in rows),
        saved_joint_frames_replayed=len(frames),
        ancestral_residue_category_pairs_checked=sum(frame['ancestral_pairs'] for frame in frames),
        native_inference_runs_started=0, existing_jobs_restarted=False,
        scientific_eligibility=False, posterior_qualified=False, source_hashes=bindings,
        scope='Every completed checkpoint present at the initial snapshot is replayed from native output; all exported array keys, shapes, dtypes and values are checked. Strict mean-one rate validation is unchanged. Partial sampler progress does not establish full completion, historical allocation repair or posterior adequacy. No original outputs are written or replaced.')
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
