#!/usr/bin/env python3
"""Check spawned output auditing, resume, and preservation of error records."""
import hashlib
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import tempfile
import numpy as np
import run_whole_protein_ml as runner
import audit_whole_protein_ml_outputs as audit
from screen_duplication_alignment_reuse import sha


def main():
    rng=np.random.default_rng(9019)
    with tempfile.TemporaryDirectory() as directory:
        root=Path(directory);factors=root/'factors';factors.mkdir()
        for i in range(5):np.savez(factors/f'tree{i}.npz',factor=rng.normal(size=(24,2))*.1)
        matrix=np.column_stack([rng.normal(size=24),np.ones(24),rng.normal(size=(24,2))]).astype('<f8')
        identities=np.asarray([hashlib.sha256(str(i).encode()).hexdigest() for i in range(24)],dtype='S64')
        spec=dict(records=24,columns=['rmsd','intercept','a','b'],values_sha256=hashlib.sha256(matrix.tobytes()).hexdigest(),ordered_identity_sha256=hashlib.sha256(identities.tobytes()).hexdigest())
        key=runner.digest(spec);path=root/'input.npz'
        np.savez(path,matrix=matrix,row_identity=identities,background=np.repeat(np.arange(6),4),family=np.repeat(np.arange(2),12),pattern_rows=np.arange(24))
        item=dict(fit_input_id=key,path=path.name,sha256=sha(path),recipe=dict(specification=spec))
        runner.initialize(factors);rows=runner.fit_case((item,str(root),str(root/'fits'),'fixture-fit-plan'))
        with ProcessPoolExecutor(max_workers=1,mp_context=multiprocessing.get_context('spawn'),initializer=audit.initialize,initargs=(str(factors),)) as pool:
            task=(item,rows,str(root),str(root/'audit'),'fixture-fit-plan','fixture-audit-plan')
            first=pool.submit(audit.audit_case,task).result()
            assert len(first)==5 and all(r['numerical_fit_verified'] for r in first)
            proof_path=root/'audit'/'inputs'/key[:2]/(key+'.json');before=proof_path.stat().st_mtime_ns
            assert pool.submit(audit.audit_case,task).result()==first and proof_path.stat().st_mtime_ns==before
            path=Path(rows[0]['path']);saved=json.loads(path.read_text())
            saved['payload']=dict(status='fit_error_requires_review',error_type='FixtureError',error='fixture',traceback='fixture traceback')
            saved['payload_sha256']=runner.digest(saved['payload']);path.write_text(json.dumps(saved))
            try:pool.submit(audit.audit_case,task).result()
            except AssertionError:pass
            else:raise AssertionError('Changed fit accepted against old binding')
            rows[0].update(sha256=sha(path),status='fit_error_requires_review')
            revised=(item,rows,str(root),str(root/'audit-error'),'fixture-fit-plan','fixture-audit-error')
            second=pool.submit(audit.audit_case,revised).result()
            assert sum(not r['numerical_fit_verified'] for r in second)==1
            assert sum(r['candidates_checked'] for r in second)==88
    result=dict(status='passed_whole_protein_output_audit_worker_fixtures',valid_tree_fits=5,candidate_likelihoods_checked=110,checks=['spawned full payload replay','unchanged resume','changed source rejected','error records explicitly unverified'],source_hashes={p:sha(p) for p in ['scripts/audit_whole_protein_ml_outputs.py','scripts/check_whole_protein_output_audit.py','scripts/check_whole_protein_fit_payload.py','scripts/run_whole_protein_ml.py']},scope='Synthetic worker and resume checks only; production output audit not launched.')
    with Path('metadata/whole_protein_output_audit_fixtures_20260929.json').open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
