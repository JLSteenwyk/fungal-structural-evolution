#!/usr/bin/env python3
"""Check full-model worker output, exact restart and corruption rejection."""
import json
import tempfile
from pathlib import Path
import numpy as np
import run_full_matched_models as runner
from screen_duplication_domain_alignment_coverage import sha


def main():
    rng=np.random.default_rng(60193);n=60
    runner.FACTORS={f'tree_{i}':rng.normal(size=(n,4))/2 for i in range(5)}
    bg=np.repeat(np.arange(15),4);family=bg//3
    matrix=rng.normal(size=(n,5));matrix[:,-1]=0
    with tempfile.TemporaryDirectory(prefix='fungal-worker-check-') as directory:
        task=('a'*64,bg,family,np.arange(n),matrix,directory,'b'*64)
        first=runner.fit_case(task)
        assert len(first)==5 and all(r['status']!='fit_error_requires_review' for r in first),first
        for row in first:
            payload=json.loads(Path(row['path']).read_text())['payload']
            assert payload['active_covariates']==[True,True,True,False]
            assert payload['direct_candidate_readback']['candidates_checked']==22
            assert np.asarray(payload['raw_unit_beta']).shape==(4,)
        assert runner.fit_case(task)==first
        path=Path(first[0]['path']);data=json.loads(path.read_text())
        data['payload']['records']+=1;path.write_text(json.dumps(data))
        try:runner.fit_case(task)
        except AssertionError:pass
        else:raise AssertionError('Corrupted saved output was accepted')
    result=dict(status='passed_full_model_worker_synthetic_restart_checks',tree_cases=5,direct_candidate_checks=110,
                resume_identical=True,corrupted_payload_rejected=True,
                pins={p:sha(p) for p in ['scripts/run_full_matched_models.py',__file__]},
                scope='Worker execution, zero-constant removal, direct candidate readback, raw coefficient shape, exact resume and corruption rejection. Computational fixtures only; full fitting-input inventory audit remains a required launch gate.')
    Path('metadata/full_matched_model_worker_checks_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
