#!/usr/bin/env python3
"""Verify independent geometry readback accepts correct examples and rejects corruption."""
import numpy as np
from assess_background_alignment_geometry import geometry
from readback_background_alignment_geometry import verify_row

rng=np.random.default_rng(776)
examples=[(np.zeros((1,3)),np.zeros((1,3))),
          (np.array([[0.,0.,0.],[1.,0.,0.]]),np.array([[0.,0.,0.],[2.,0.,0.]])),
          (np.eye(3),np.eye(3)),
          (np.vstack([np.eye(3),-np.eye(3)]),np.vstack([np.eye(3),-np.eye(3)])*[-1,1,1]),
          (rng.normal(size=(60,3)),rng.normal(size=(60,3)))]
for x,y in examples:
    row={k:str(v) for k,v in geometry(x,y).items()}
    verify_row(row,x,y)
    for column,value in [('minimum_rotation_curvature','10000000'),('rank_left','9'),('geometry_status','invalid'),('relative_rotation_curvature','-999')]:
        corrupt=dict(row);corrupt[column]=value
        try:verify_row(corrupt,x,y)
        except ValueError:pass
        else:raise AssertionError('Corruption accepted: '+column)
print('Quaternion readback passed point, two-point, planar, reflected and random cases; altered curvature, rank, status and ratio rejected.')
