"""Construct a nonnegative support certificate; do not claim an LP optimum."""
import math
import numpy as np


def project_weights(x, saved):
    x=np.asarray(x,dtype=float)
    assert x.ndim==2 and len(x)==saved['records'] and np.isfinite(x).all()
    scales=np.max(abs(x),axis=0);scales[scales==0]=1.
    np.testing.assert_array_equal(scales,saved['scales'])
    indices=saved.get('support_indices',[]);weights=saved.get('support_weights',[])
    positive=[(i,float(w)) for i,w in zip(indices,weights) if w>0]
    if not positive:return dict(classification='unresolved_no_positive_weights')
    total=math.fsum(w for _,w in positive)
    indices=[i for i,_ in positive];weights=[w/total for _,w in positive]
    z=x/scales
    barycenter=np.asarray([math.fsum(w*float(z[i,j]) for i,w in zip(indices,weights)) for j in range(x.shape[1])])
    np.testing.assert_allclose(np.asarray(weights)@z[indices],barycenter,rtol=1e-10,atol=1e-12)
    distance=float(max(abs(barycenter)))
    assert all(w>=0 for w in weights) and abs(math.fsum(weights)-1)<=1e-12
    return dict(classification='supported_by_projected_nonnegative_weights' if distance<=1e-8 else 'unresolved_after_projection',
        records=len(x),scales=scales.tolist(),support_indices=indices,support_weights=weights,
        barycenter=barycenter.tolist(),primal_distance=distance,
        removed_negative_mass=math.fsum(-w for w in saved.get('support_weights',[]) if w<0),
        method='Set negative saved weights to zero, normalize positive weights, recompute against original matrix. No optimizer or optimum claim.')
