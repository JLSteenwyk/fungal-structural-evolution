#!/usr/bin/env python3
"""Exercise real spawned workers, immutable resume and corruption rejection."""
import hashlib
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import tempfile
import numpy as np
from run_whole_protein_ml import initialize, fit_case, digest
from screen_duplication_alignment_reuse import sha


def main():
    rng = np.random.default_rng(9131)
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        factors = root/'factors'
        factors.mkdir()
        for i in range(5):
            np.savez(factors/f'tree{i}.npz', factor=rng.normal(size=(24, 2))*.1)
        values = np.column_stack([rng.normal(size=24), np.ones(24), rng.normal(size=(24, 2))]).astype('<f8')
        identities = np.array([hashlib.sha256(str(i).encode()).hexdigest() for i in range(24)], dtype='S64')
        spec = dict(records=24, columns=['rmsd', 'intercept', 'a', 'b'], values_sha256=hashlib.sha256(values.tobytes()).hexdigest(), ordered_identity_sha256=hashlib.sha256(identities.tobytes()).hexdigest())
        identifier = digest(spec)
        path = root/'input.npz'
        np.savez(path, matrix=values, row_identity=identities, background=np.repeat(np.arange(6),4), family=np.repeat(np.arange(2),12), pattern_rows=np.arange(24))
        item = dict(fit_input_id=identifier, recipe=dict(specification=spec), path=path.name, sha256=sha(path))
        task = item, str(root), str(root/'output'), 'fixture-plan'
        with ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context('spawn'), initializer=initialize, initargs=(str(factors),)) as pool:
            first = pool.submit(fit_case, task).result()
            assert len(first)==5 and all(r['status']!='fit_error_requires_review' for r in first)
            modified = {r['path']:Path(r['path']).stat().st_mtime_ns for r in first}
            second = pool.submit(fit_case, task).result()
            assert first == second
            assert modified == {p:Path(p).stat().st_mtime_ns for p in modified}
            bad = dict(item, sha256='0'*64)
            try:
                pool.submit(fit_case, (bad, *task[1:])).result()
            except AssertionError:
                pass
            else:
                raise AssertionError('Changed source binding accepted')
            result_path = Path(first[0]['path'])
            saved = json.loads(result_path.read_text())
            saved['payload']['ratios'][0] += 1
            result_path.write_text(json.dumps(saved))
            try:
                pool.submit(fit_case, task).result()
            except AssertionError:
                pass
            else:
                raise AssertionError('Corrupt resumed payload accepted')
    result = dict(status='passed_whole_protein_ml_runner_worker_fixtures', tree_fits=5,
                  checks=['spawned worker execution', 'five tree outputs', 'byte-identical and mtime-preserving resume', 'source hash mismatch rejected', 'resumed payload corruption rejected'],
                  source_hashes={p:sha(p) for p in ['scripts/run_whole_protein_ml.py', 'scripts/fit_whole_protein_ml.py', 'scripts/check_whole_protein_ml_runner.py']},
                  scope='Worker and restart behavior on synthetic inputs only. Full-data gates and allocation remain pending; no production run launched.')
    with Path('metadata/whole_protein_ml_runner_fixtures_20260929.json').open('x') as handle:
        json.dump(result,handle,indent=2)
        handle.write('\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    main()
