#!/usr/bin/env python3
"""Known quadratic optima and rejection paths for lower-face proposals."""
import json
from pathlib import Path
import numpy as np
from lower_face_gradient_polish import propose
from screen_duplication_alignment_reuse import sha


def main():
    passed=[]
    for optimum in (np.array([.2,-1.,.4]), np.array([-1.,-1.,.4])):
        initial=np.array([.3,0.,.5]) if optimum[0]>0 else np.array([0.,0.,.5])
        objective=lambda x:float(np.sum((x-optimum)**2))
        gradient=lambda x:2*(x-optimum)
        result=propose(objective,gradient,initial,2.)
        assert result['status']=='two_lower_face_proposals_passed'
        for trial in result['proposals']:
            np.testing.assert_allclose(trial['theta'],np.maximum(optimum,0),atol=1e-10)
        passed.append('known_face_optimum_'+str(len(result['fixed_axes'])))
    obj=lambda x:float((x[0]-.2)**2+(x[1]+1)**2+(x[2]-.4)**2)
    grad=lambda x:2*(x-np.array([.2,-1.,.4]))
    cases=[([.3,0,.5],lambda x:np.array([1.,-1.,1.]),'boundary_sign_requires_review'),([0,0,0],grad,'not_applicable_without_mixed_lower_face'),([.3,.2,.5],grad,'not_applicable_without_mixed_lower_face'),([1e-7,0,.5],grad,'not_applicable_near_other_boundary'),([.2,0,.4],grad,'no_free_gradient_correction_needed')]
    for t,g,status in cases:
        assert propose(obj,g,t,2.)['status']==status
        passed.append(status)
    for invalid in ([np.nan,0,.5],[-1,0,.5],[3,0,.5],[1,2]):
        try:propose(obj,grad,invalid,2.)
        except ValueError:pass
        else:raise AssertionError('invalid parameters accepted')
    bad=propose(obj,lambda x:grad(x)+np.array([.05,0,0]),[.3,0,.5],2.)
    assert bad['status']=='lower_face_proposals_require_review'
    assert not any(t['checks']['direct_free_gradients_pass'] for t in bad['proposals'])
    passed.append('inconsistent_gradient_rejected_by_direct_differences')
    proof=dict(status='passed_lower_face_polish_fixtures',cases=passed,invalid_parameter_cases=4,helper_sha256=sha('scripts/lower_face_gradient_polish.py'),checker_sha256=sha(__file__))
    with Path('metadata/lower_face_polish_fixtures_20260929.json').open('x') as h:
        json.dump(proof,h,indent=2);h.write('\n')
    print(json.dumps(proof))


if __name__=='__main__':
    main()
