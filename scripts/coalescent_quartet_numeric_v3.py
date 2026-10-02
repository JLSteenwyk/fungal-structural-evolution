"""Exact native ASTRAL 5.7.8 effective-N rule with strict numerical checks.

The native code begins with the number of genes with all four incident clades
available. It substitutes the fractional resolved evidence only when their
absolute difference exceeds 0.001. This is an algorithmic rule, not a relaxed
comparison tolerance. Raw fractional evidence remains separately recorded.
Reference: upstream WQInference.java scoreBranches; actual pinned-jar boundary
contracts are required because upstream master alone does not prove behavior.
"""
import math
import numpy as np
from scipy.special import betaln,betaincc,logsumexp
from coalescent_quartet_audit_v2 import close,map_length


def effective_n(counts,available_genes):
    assert isinstance(available_genes,(int,np.integer)) and available_genes>=0
    counts=np.asarray(counts,dtype=float)
    assert counts.shape==(3,) and np.isfinite(counts).all() and min(counts)>=0
    total=float(counts.sum());assert total<=available_genes+1e-8
    return total if abs(total-available_genes)>.001 else float(available_genes)


def posterior(counts,total,prior=.5):
    counts=np.asarray(counts,dtype=float)
    assert np.isfinite(counts).all() and min(counts)>=0 and total>=0 and prior>0
    assert max(counts)<=total+1e-8
    a,b=counts+1,total-counts+2*prior
    logs=counts*math.log(2)+betaln(a,b)+np.log(betaincc(a,b,1/3))
    result=np.exp(logs-logsumexp(logs));assert np.isfinite(result).all()
    return result


def numerical_check(values,length,independent_counts,available_genes,tolerance=2e-8):
    native=np.array([values['f1'],values['f2'],values['f3']])
    independent=np.asarray(independent_counts,dtype=float)
    assert np.isfinite(native).all() and min(native)>=0
    close(native[0],independent[0],tolerance)
    for a,b in zip(sorted(native[1:]),sorted(independent[1:])):close(a,b,tolerance)
    # Use checked independent evidence, oriented to the native alternative labels.
    alternatives=sorted(independent[1:],reverse=native[1]>native[2])
    aligned=np.array([independent[0],*alternatives])
    total=float(independent.sum());en=effective_n(aligned,available_genes)
    close(values['EN'],en,tolerance)
    expected_pp=posterior(aligned,en)
    for i in range(3):
        if en==0:assert math.isnan(values['q'+str(i+1)])
        else:close(values['q'+str(i+1)],aligned[i]/en,tolerance)
        close(values['pp'+str(i+1)],expected_pp[i],tolerance)
    expected_length=map_length(independent[0],en)
    close(length,expected_length,tolerance)
    close(sum(values['pp'+str(i+1)] for i in range(3)),1,tolerance)
    return dict(independent_effective_genes=total,independent_available_genes=int(available_genes),
        independent_native_effective_genes=en,native_effective_n_adjustment_applied=abs(total-available_genes)>.001,
        native_and_independent_counts=native.tolist(),independent_counts_canonical_groups=independent.tolist(),
        independent_posterior=expected_pp.tolist(),independent_map_length=expected_length,q_available=en>0)
