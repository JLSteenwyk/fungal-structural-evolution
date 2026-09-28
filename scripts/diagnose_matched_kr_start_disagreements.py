"""Recompute full-face candidate scores for every unresolved synthetic refit."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha,write_json
from cached_matched_likelihood import CachedMatchedLikelihood
from matched_reml_gradient import evaluate_gradient
from matched_mixed_covariance import MatchedCovariance,profiled_reml
from simulate_matched_working_model import simulate


def main():
    threadpool_limits(1)
    root=Path('results/model_validation/matched-kr-simulations-20260928-v1')
    receipt=json.loads((root/'receipt.json').read_text())
    assert sha(root/'manifest.json')==receipt['artifacts']['manifest.json']
    manifest=json.loads((root/'manifest.json').read_text())
    selected=[r for r in manifest if r['refit_status']=='refit_requires_review']
    assert len(selected)==184
    rows=[];bindings={};maximum_error=0.;arrays={}
    for entry in selected:
        path=Path(entry['path']);assert sha(path)==entry['sha256'];bindings[str(path)]=sha(path)
        saved=json.loads(path.read_text());case=saved['case'];r=saved['refit'];f=r['fit']
        if case['design'] not in arrays:
            with np.load(case['design'],allow_pickle=False) as h:arrays[case['design']]={k:h[k] for k in h.files}
        a=arrays[case['design']]
        y=simulate(a['background'],a['family'],a['factor'],a['design'],a['beta'],case['scale'],case['ratios'],1,
                   np.random.default_rng(np.random.SeedSequence(r['seed_entropy'])))[:,0]
        assert hashlib.sha256(np.asarray(y,dtype='<f8').tobytes()).hexdigest()==r['response_sha256']
        cache=CachedMatchedLikelihood(a['background'],a['family'],a['factor'],a['design'],y)
        upper=np.log1p(f['maximum_ratio'])
        for candidate in f['candidates']:
            if not all(candidate['active']) or candidate['start']=='original_unoptimized_reference':continue
            theta=np.array(candidate['theta']);ratios=np.expm1(theta)
            value=evaluate_gradient(cache,ratios);g=value['log1p_ratio_gradient']
            pg=np.where(theta<=1e-7,np.minimum(g,0),np.where(theta>=upper-1e-7,np.maximum(g,0),g))
            direct=profiled_reml(MatchedCovariance(a['background'],a['family'],a['factor'],1.,*ratios),a['design'],y)
            np.testing.assert_allclose(direct['negative_profiled_reml'],candidate['objective'],rtol=1e-9,atol=1e-7)
            np.testing.assert_allclose(value['negative_profiled_reml'],candidate['objective'],rtol=1e-9,atol=1e-7)
            maximum_error=max(maximum_error,abs(direct['negative_profiled_reml']-candidate['objective']))
            gap=candidate['objective']-f['negative_profiled_reml'];norm=float(np.max(abs(pg)))
            diagnosis=('agrees_with_selected' if gap<=1e-5 else
                       'higher_objective_nonstationary' if norm>1e-3 else 'higher_objective_stationary_requires_review')
            rows.append(dict(case=case['id'],replicate=saved['replicate'],start=candidate['start'],
                             original_success=candidate['success'],original_message=candidate['message'],
                             objective_gap=gap,projected_gradient_maximum=norm,diagnosis=diagnosis,
                             source_path=str(path),source_sha256=entry['sha256']))
    assert len(rows)==184*4
    frame=pd.DataFrame(rows);assert not frame.duplicated(['case','replicate','start']).any()
    out=Path('results/model_validation/matched-kr-start-diagnostics-20260928-v1');out.mkdir(parents=True,exist_ok=False)
    frame.to_csv(out/'full_face_candidates.tsv',sep='\t',index=False)
    counts=dict(Counter(frame.diagnosis))
    result=dict(status='completed_all_184_review_fit_full_face_score_diagnostics',fits=184,candidates=len(rows),
        diagnosis_counts=counts,maximum_direct_objective_disagreement=maximum_error,
        nonstationary_higher_candidates_reported_success=int(((frame.diagnosis=='higher_objective_nonstationary')&frame.original_success).sum()),
        source_receipt_sha256=sha(root/'receipt.json'),source_files=bindings,
        pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/matched_reml_gradient.py'),
            Path('scripts/cached_matched_likelihood.py'),Path('scripts/matched_mixed_covariance.py'),
            Path('scripts/simulate_matched_working_model.py')]},
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='All four full-face starts in every unresolved synthetic refit. Direct likelihood '
              'and analytic projected-score replay only; no optimizer restart or qualification '
              'change. Stationarity is a local diagnostic, not a global-optimum proof.')
    write_json(out/'receipt.json',result);write_json(Path('metadata/matched_kr_start_diagnostics_20260928.json'),result)
    print({k:v for k,v in result.items() if k not in ['source_files','pins','artifacts']})


if __name__=='__main__':main()
