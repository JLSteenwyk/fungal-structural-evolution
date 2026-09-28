#!/usr/bin/env python3
"""Lossless exact trace grouping; pattern equality is not biological homology."""
import numpy as np


def group_patterns(values):
    """Input: chain, saved iteration, then one or more coordinate axes.

Return distinct full temporal traces and a pattern id for every coordinate.
No sorting of draws, relabeling of states, thinning or frequency binning.
"""
    values=np.asarray(values)
    if values.dtype!=np.uint8 or values.ndim<3 or any(n==0 for n in values.shape):
        raise ValueError('Nonempty uint8 chain/draw/coordinate array required')
    shape=values.shape
    flat=values.reshape(shape[0],shape[1],-1)
    mapping=np.empty(flat.shape[2],dtype=np.int32)
    lookup={};keys=[]
    for coordinate in range(flat.shape[2]):
        key=flat[:,:,coordinate].tobytes(order='C')
        index=lookup.get(key)
        if index is None:
            index=len(keys);lookup[key]=index;keys.append(key)
        mapping[coordinate]=index
    patterns=np.frombuffer(b''.join(keys),dtype=np.uint8).reshape(len(keys),shape[0],shape[1]).copy()
    return patterns,mapping.reshape(shape[2:])


def reconstruct_patterns(patterns,mapping):
    patterns=np.asarray(patterns);mapping=np.asarray(mapping)
    if patterns.ndim!=3 or patterns.dtype!=np.uint8 or mapping.ndim<1 or not np.issubdtype(mapping.dtype,np.integer):
        raise ValueError('Invalid patterns or coordinate map')
    if mapping.size==0 or np.any(mapping<0) or np.any(mapping>=len(patterns)):
        raise ValueError('Invalid pattern ids')
    return np.moveaxis(patterns[mapping],(-2,-1),(0,1))
