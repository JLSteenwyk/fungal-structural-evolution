"""Reproduce the exact source-bound 60-fit replay fixture set."""
import json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha
from audit_matched_simulation_cache import replay


def main():
    threadpool_limits(1)
    receipt=json.loads(Path('metadata/matched_simulation_cache_replay_checks_20260928.json').read_text())
    assert sha('scripts/audit_matched_simulation_cache.py')==receipt['script_sha256']
    inputs={};fits=[]
    for name,digest in receipt['source_bindings'].items():
        path=Path(name);assert sha(path)==digest
        if path.suffix=='.npz':
            with np.load(path,allow_pickle=False) as h:inputs[path.stem]={k:h[k] for k in h.files}
        else:fits.append(json.loads(path.read_text()))
    cache=Path('results/model_validation/matched-simulation-input-cache-20260928-v1')
    production=json.loads(Path('metadata/full_matched_working_model_plan_20260927.json').read_text())
    factors={}
    for path in (cache/'factors').glob('*.npz'):
        source=str(Path(production['factors'])/path.name)
        assert sha(path)==production['pins'][source]
        with np.load(path,allow_pickle=False) as h:factors[path.stem]=h['factor']
    errors=[]
    for saved in fits:
        errors.append(replay(inputs[saved['fit_input_id']],saved['payload'],factors[saved['tree']])['objective_absolute_error'])
    assert len(inputs)==12 and len(fits)==60
    np.testing.assert_allclose(max(errors),receipt['maximum_objective_error'],rtol=1e-5,atol=1e-11)
    bad={k:v.copy() for k,v in inputs[saved['fit_input_id']].items()};bad['matrix'][0,0]+=1
    try:replay(bad,saved['payload'],factors[saved['tree']])
    except AssertionError:pass
    else:raise AssertionError('Changed response accepted')
    print('Reproduced all60 source-bound fit checks; changed response rejected')


if __name__=='__main__':main()
