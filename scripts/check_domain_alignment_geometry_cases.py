#!/usr/bin/env python3
"""Check rigid-fit degeneracy and curvature against independent finite differences."""
import numpy as np
from assess_domain_alignment_geometry import geometry

line=np.array([[0.,0.,0.],[1.,0.,0.],[2.,0.,0.]])
plane=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
for x in [np.zeros((1,3)),line[:2],line]:
    assert geometry(x,x)['geometry_status']=='degenerate_at_numeric_tolerance'
assert geometry(plane,plane)['geometry_status']=='unique_at_numeric_tolerance'
# Isotropic reflected correspondences have no unique best proper rotation.
x=np.vstack([np.eye(3),-np.eye(3)]);y=x*np.array([-1,1,1])
assert geometry(x,y)['geometry_status']=='degenerate_at_numeric_tolerance'
rng=np.random.default_rng(12871);x=rng.normal(size=(40,3));y=x@np.diag([1.,2.,-3.])+rng.normal(size=(40,3))*.02
x-=x.mean(0);y-=y.mean(0)
u,_,vt=np.linalg.svd(x.T@y);d=np.diag([1.,1.,np.linalg.det(u@vt)]);rotation=u@d@vt

def exp_rotation(v):
    angle=np.linalg.norm(v)
    if angle==0:return np.eye(3)
    a=v/angle;k=np.array([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0]])
    return np.eye(3)+np.sin(angle)*k+(1-np.cos(angle))*(k@k)

def cost(v):return np.sum((x@rotation@exp_rotation(v)-y)**2)

step=1e-3;axes=np.eye(3)*step;hessian=np.zeros((3,3));base=cost(np.zeros(3))
for i in range(3):
    hessian[i,i]=(cost(axes[i])+cost(-axes[i])-2*base)/step**2
    for j in range(i):
        hessian[i,j]=hessian[j,i]=(cost(axes[i]+axes[j])-cost(axes[i]-axes[j])-cost(-axes[i]+axes[j])+cost(-axes[i]-axes[j]))/(4*step**2)
g=geometry(x,y)
assert abs(np.linalg.eigvalsh(hessian)[0]/2-g['minimum_rotation_curvature'])<2e-4
h=geometry(y,x)
for key in ['cross_s1','cross_s2','cross_s3','minimum_rotation_curvature','relative_rotation_curvature']:
    assert np.isclose(g[key],h[key],rtol=1e-12,atol=1e-12)
# Independent rotation and translation preserve spectra and curvature.
q=exp_rotation(np.array([.4,-.2,.7]));t=geometry(x@q+[100,-30,40],y+[-40,2,80])
for key in ['cross_s1','cross_s2','cross_s3','minimum_rotation_curvature']:
    assert np.isclose(g[key],t[key],rtol=1e-12,atol=1e-12)
print('Point, line, plane, reflection degeneracy, finite-difference curvature, order symmetry and rigid-motion invariance passed.')
