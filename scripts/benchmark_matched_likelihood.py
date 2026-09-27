#!/usr/bin/env python3
"""Numerical runtime measurements at observed full-design sizes, not data fits."""
import json
import os
import platform
import time
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_info, threadpool_limits
from matched_mixed_covariance import MatchedCovariance, profiled_reml
from screen_duplication_domain_alignment_coverage import sha


def main():
    source=Path('results/structural_comparisons/matched-domain-record-summaries-20260927-v1/record_summary.tsv')
    proof=Path('metadata/matched_domain_record_summary_completed_20260927.json')
    receipt=Path('results/structural_comparisons/matched-domain-record-summaries-20260927-v1/receipt.json')
    assert sha(source)==json.loads(receipt.read_text())['artifacts'][source.name]
    assert json.loads(proof.read_text())['source_receipt_sha256']==sha(receipt)
    d=pd.read_csv(source,sep='\t')
    median=float(d.matched_records.median())
    sizes=[int(d.matched_records.min()),int(d.loc[d.matched_records<=median,'matched_records'].max()),
           int(d.loc[d.matched_records>=median,'matched_records'].min()),int(d.matched_records.max())]
    labels=['minimum','median_lower_observed','median_upper_observed','maximum']
    rank=242
    rng=np.random.default_rng(442019)
    measurements=[]
    with threadpool_limits(limits=1):
        pools=threadpool_info()
        assert all(p['num_threads']==1 for p in pools)
        for label,size in zip(labels,sizes):
            row=d.iloc[(d.matched_records-size).abs().argmin()]
            n,g,f=int(row.matched_records),int(row.backgrounds),int(row.families)
            bg=np.concatenate([np.arange(g),rng.integers(0,g,n-g)])
            rng.shuffle(bg)
            bg_family=np.concatenate([np.arange(f),rng.integers(0,f,g-f)])
            rng.shuffle(bg_family)
            family=bg_family[bg]
            factor=rng.normal(size=(n,rank))/np.sqrt(rank)
            x=np.column_stack([np.ones(n),rng.normal(size=(n,4))])
            y=rng.normal(size=n)
            # Warm up imports and kernels outside timing, then repeat fresh construction.
            profiled_reml(MatchedCovariance(bg,family,factor,1.,.2,.5,.3),x,y)
            for ratios in [(0.,0.,0.),(.2,.5,0.),(.2,.5,.3)]:
                for repeat in range(5):
                    cpu=time.process_time();wall=time.perf_counter()
                    covariance=MatchedCovariance(bg,family,factor,1.,*ratios)
                    value=profiled_reml(covariance,x,y)
                    elapsed=time.perf_counter()-wall;used=time.process_time()-cpu
                    assert np.isfinite(value['negative_profiled_reml'])
                    measurements.append(dict(size_label=label,records=n,backgrounds=g,families=f,
                                             factor_columns=rank,fixed_columns=5,ratios=list(ratios),
                                             repeat=repeat,cpu_seconds=used,wall_seconds=elapsed))
            print('Measured numerical likelihood size',label,n,flush=True)
    optimizer_root=Path('results/model_validation/matched-mixed-optimizer-20260927-v1')
    optimizer_receipt=json.loads((optimizer_root/'receipt.json').read_text())
    for name,digest in optimizer_receipt['artifacts'].items():assert sha(optimizer_root/name)==digest
    attempts=[]
    for label in ['residual_only','mixed','strong_background']:
        fit=json.loads((optimizer_root/(label+'.json')).read_text())['fit']
        attempts.append(dict(case=label,optimizer_objective_calls=sum(c['evaluations'] for c in fit['candidates']),
                             additional_readback_and_gradient_calls=7))
    summary=[]
    for label in labels:
        for ratios in [(0.,0.,0.),(.2,.5,0.),(.2,.5,.3)]:
            selected=[m for m in measurements if m['size_label']==label and m['ratios']==list(ratios)]
            summary.append(dict(size_label=label,records=selected[0]['records'],ratios=list(ratios),
                                median_cpu_seconds=float(np.median([m['cpu_seconds'] for m in selected])),
                                median_wall_seconds=float(np.median([m['wall_seconds'] for m in selected])),
                                minimum_wall_seconds=min(m['wall_seconds'] for m in selected),
                                maximum_wall_seconds=max(m['wall_seconds'] for m in selected)))
    assert len(measurements)==60 and len(summary)==12
    result=dict(status='complete_synthetic_likelihood_runtime_measurements',measurements=measurements,summaries=summary,
                observed_record_count_median=median,benchmark_record_counts=sizes,
                observed_synthetic_optimizer_calls=attempts,host=platform.node(),platform=platform.platform(),
                numpy_version=np.__version__,threadpools=pools,logical_cpus=os.cpu_count(),
                pins={str(p):sha(p) for p in [source,proof,receipt,Path(__file__),Path('scripts/matched_mixed_covariance.py'),optimizer_root/'receipt.json']+[optimizer_root/(label+'.json') for label in ['residual_only','mixed','strong_background']]},
                scope='Synthetic numerical inputs at observed minimum, lower/upper observations bracketing the median, and maximum record/group counts; full 242-column factor and five fixed columns. Five repeats at each of three component settings. No biological subset fit, optimizer convergence forecast, confidence calibration or whole-grid ETA; real rank and conditioning may differ, and observed wall times reflect competing live jobs.')
    Path('metadata/matched_likelihood_runtime_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(summaries=summary,optimizer_calls=attempts),indent=2))


if __name__=='__main__':main()
