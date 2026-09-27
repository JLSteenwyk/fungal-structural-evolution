#!/usr/bin/env python3
"""Exhaustively verify partition sums, marginals and MAP for small interval sets."""
import itertools,json
from pathlib import Path
import numpy as np
from compatible_gap_distribution import CompatibleGaps
from prepare_case_ancestral_neighborhoods import sha


def valid(intervals,bits):
    selected=[i for i,x in enumerate(bits) if x]
    return all(max(intervals[i][0],intervals[j][0])>min(intervals[i][1],intervals[j][1])+1 for i,j in itertools.combinations(selected,2))


def main():
    rng=np.random.default_rng(20260927);tests=[]
    fixed=[([(1,3),(2,4),(4,6),(8,9)],[.4,.4,.4,.7]), ([(1,3),(2,4)],[1.,1.]), ([(1,3),(8,9)],[1.,0.]),([],[])]
    for _ in range(80):
        starts=rng.choice(30,7,replace=False);intervals=[(int(a),int(a+rng.integers(0,9))) for a in starts]
        fixed.append((intervals,rng.uniform(.001,.999,7)))
    for intervals,p in fixed:
        p=np.array(p);model=CompatibleGaps(intervals,p);weights=[];states=[]
        for bits in itertools.product([0,1],repeat=len(p)):
            if valid(intervals,bits):
                states.append(bits);weights.append(float(np.prod(np.where(bits,p,1-p))))
        total=sum(weights)
        assert model.feasible==(total>0)
        if total:
            expected=np.array(weights)@np.array(states).reshape(len(states),len(p))/total
            got=model.marginals();error=float(np.max(abs(got-expected))) if len(p) else 0.
            assert error<1e-12 and abs(model.log_compatibility_probability-np.log(total))<1e-12
            chosen=model.configuration();bits=[int(i in chosen) for i in range(len(p))]
            assert valid(intervals,bits)
            w=float(np.prod(np.where(bits,p,1-p)));assert abs(w-max(weights))<1e-12
            for _ in range(50):
                draw=model.configuration(rng);assert valid(intervals,[int(i in draw) for i in range(len(p))])
        else:
            error=None
            try:model.marginals()
            except ValueError:pass
            else:raise AssertionError('infeasible distribution accepted')
        tests.append(dict(intervals=len(p),compatible_mass=total,maximum_marginal_error=error))
    result=dict(status='all_exhaustive_compatible_gap_distribution_checks_passed',cases=tests,
        pins={p:sha(p) for p in [__file__,'scripts/compatible_gap_distribution.py']},
        scope='Exact DP for product-Bernoulli weights conditioned on nonoverlap/nonadjacency, including zero/one weights and infeasible support. This is a declared surrogate that changes marginals, not a fitted joint indel evolutionary model. Random draws checked for compatibility; sampling frequencies not validated by these tests.')
    Path('metadata/compatible_gap_distribution_validation_20260927.json').write_text(json.dumps(result,indent=2)+'\n');print(result['status'],len(tests))

if __name__=='__main__':main()
