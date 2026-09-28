#!/usr/bin/env python3
"""Bounded synthetic timing/serialization benchmark, not biological inference."""
import argparse
import gzip
import json
from pathlib import Path
import time
import numpy as np
from ancestral_chain_diagnostics import sha
from ancestral_categorical_diagnostics import diagnose_states


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();assert not args.output.exists()
    rng=np.random.default_rng(20260928);alphabet='ACDEFGHIKLMNPQRSTVWYX-'
    # Exclude import/first-call overhead from per-pattern arithmetic scaling.
    diagnose_states(rng.integers(0,2,size=(4,75)),alphabet)
    rows=[];started=time.monotonic()
    for draws in [50,75]:
        for kind in ['constant','rare','sticky2','iid2','iid4','iid22']:
            times=[];payload=[]
            for repetition in range(10):
                if kind=='constant':values=np.zeros((4,draws),dtype=np.uint8)
                elif kind=='rare':
                    values=np.zeros((4,draws),dtype=np.uint8);values[0,repetition]=1
                elif kind=='sticky2':values=np.repeat(rng.integers(0,2,size=(4,(draws+9)//10)),10,axis=1)[:,:draws]
                else:values=rng.integers(0,int(kind[3:]),size=(4,draws))
                before=time.monotonic();result=diagnose_states(values,alphabet);times.append(time.monotonic()-before)
                payload.append(json.dumps(dict(pattern_id=repetition,coordinate_multiplicity=1,diagnostic=result),allow_nan=False,separators=(',',':'))+'\n')
            raw=''.join(payload).encode();compressed=gzip.compress(raw,mtime=0)
            rows.append(dict(draws_per_chain=draws,fixture=kind,repetitions=10,seconds_per_pattern=times,median_seconds=float(np.median(times)),maximum_seconds=max(times),raw_bytes_per_pattern=len(raw)/10,gzip_bytes_per_pattern=len(compressed)/10))
    files=[__file__,'scripts/ancestral_categorical_diagnostics.py','scripts/ancestral_chain_diagnostics.py','environments/ancestral-diagnostics-20260927.lock.txt']
    result=dict(status='bounded_synthetic_runtime_benchmark_complete',patterns=120,elapsed_seconds=time.monotonic()-started,rows=rows,pins={p:sha(p) for p in files},scope='Synthetic performance fixtures only. Real unique-pattern counts, compression and heterogeneity differ; no convergence or full-grid ETA claim.')
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['patterns','elapsed_seconds']}));print('Max measured seconds per pattern:',max(r['maximum_seconds'] for r in rows))


if __name__=='__main__':main()
