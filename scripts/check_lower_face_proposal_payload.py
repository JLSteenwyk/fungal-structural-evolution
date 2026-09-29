#!/usr/bin/env python3
"""Replay lower-face payload without calling its proposal implementation."""
import numpy as np


def check_proposal(result, objective, gradient, theta, upper):
    theta=np.asarray(theta,dtype=float)
    fixed=np.where(theta==0)[0];free=np.where(theta>0)[0]
    if len(fixed)==0 or len(free)==0:
        assert result==dict(status='not_applicable_without_mixed_lower_face',proposals=[])
        return
    if min(theta[free])<=1e-5 or max(theta[free])>=upper-1e-5:
        assert result==dict(status='not_applicable_near_other_boundary',proposals=[])
        return
    base=objective(theta);g=np.asarray(gradient(theta))
    np.testing.assert_allclose(result['original_gradient'],g,rtol=1e-9,atol=1e-9)
    if min(g[fixed])<0:
        assert result['status']=='boundary_sign_requires_review' and result['proposals']==[]
        return
    if max(abs(g[free]))<=1e-3:
        assert result['status']=='no_free_gradient_correction_needed' and result['proposals']==[]
        return
    np.testing.assert_array_equal(result['original_theta'],theta)
    np.testing.assert_allclose(result['original_objective'],base,rtol=1e-12,atol=1e-9)
    assert result['fixed_axes']==fixed.tolist() and result['free_axes']==free.tolist() and result['upper']==upper
    assert len(result['proposals'])==2
    passed=[]
    for step,trial in zip((1e-5,1e-6),result['proposals']):
        assert trial['hessian_step']==step
        raw=np.zeros((len(free),len(free)))
        for col,j in enumerate(free):
            high=theta.copy();low=theta.copy();high[j]+=step;low[j]-=step
            raw[:,col]=(np.asarray(gradient(high))[free]-np.asarray(gradient(low))[free])/(2*step)
        hessian=(raw+raw.T)*.5;eigenvalues=np.linalg.eigvalsh(hessian)
        np.testing.assert_allclose(trial['hessian'],hessian,rtol=1e-9,atol=1e-8)
        np.testing.assert_allclose(trial['eigenvalues'],eigenvalues,rtol=1e-9,atol=1e-8)
        if not np.isfinite(eigenvalues).all() or min(eigenvalues)<=0:
            assert not trial['passed'] and trial['reason']=='nonpositive_or_nonfinite_curvature'
            passed.append(False);continue
        point=theta.copy();point[free]=theta[free]-np.linalg.solve(hessian,g[free])
        np.testing.assert_allclose(trial['theta'],point,rtol=1e-12,atol=1e-12)
        if not np.isfinite(point).all() or min(point[free])<=1e-6 or max(point[free])>=upper-1e-6:
            assert not trial['passed'] and trial['reason']=='proposal_not_interior_on_free_face'
            passed.append(False);continue
        value=objective(point);derivative=np.asarray(gradient(point))
        np.testing.assert_allclose(trial['objective'],value,rtol=1e-12,atol=1e-9)
        np.testing.assert_allclose(trial['gradient'],derivative,rtol=1e-9,atol=1e-9)
        assert len(trial['finite_differences'])==2
        free_pass=[];boundary_pass=[]
        for h,record in zip((1e-6,1e-7),trial['finite_differences']):
            assert record['step']==h
            estimates=[]
            for j in range(3):
                plus=point.copy();plus[j]+=h
                if j in fixed:
                    twice=point.copy();twice[j]+=2*h
                    estimates.append((-3*value+4*objective(plus)-objective(twice))/(2*h))
                else:
                    minus=point.copy();minus[j]-=h
                    estimates.append((objective(plus)-objective(minus))/(2*h))
            fd=np.asarray(estimates)
            np.testing.assert_allclose(record['gradient'],fd,rtol=1e-9,atol=1e-9)
            fp=bool(max(abs(fd[free]))<=1e-3);bp=bool(min(fd[fixed])>=0)
            assert record['free_gradient_pass']==fp and record['boundary_sign_pass']==bp
            free_pass.append(fp);boundary_pass.append(bp)
        checks=dict(objective_not_worsened=bool(value<=base),analytic_free_gradient_pass=bool(max(abs(derivative[free]))<=1e-3),analytic_boundary_sign_pass=bool(min(derivative[fixed])>=0),direct_free_gradients_pass=all(free_pass),direct_boundary_signs_pass=all(boundary_pass))
        assert trial['checks']==checks and trial['passed']==all(checks.values())
        passed.append(all(checks.values()))
    assert result['status']==('two_lower_face_proposals_passed' if all(passed) else 'lower_face_proposals_require_review')
