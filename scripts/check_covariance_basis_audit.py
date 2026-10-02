#!/usr/bin/env python3
"""Independent dense error-contrast tests of the sparse covariance basis audit."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from unittest.mock import patch
import numpy as np
from scipy import sparse
from scipy.linalg import null_space
from covariance_basis_audit import audit_covariance_bases


def dense_check(labels, incidence, factor, diagonal, x):
    audit = audit_covariance_bases(labels, incidence, factor, diagonal, x)
    kernels = [np.diag(diagonal)]
    kernels.extend(np.asarray(z @ z.T) if not sparse.issparse(z) else (z @ z.T).toarray()
                   for z in incidence.values())
    kernels.append(factor @ factor.T)
    # Balance the dense oracle independently in extended precision: extreme
    # units must not cause null_space to discard a scientifically active term.
    balanced = np.asarray(x, dtype=np.longdouble)
    balanced /= np.max(abs(balanced), axis=0)
    balanced /= np.sqrt(np.sum(balanced*balanced, axis=0))
    contrasts = null_space(np.asarray(balanced.T, dtype=float))
    restricted = [contrasts.T @ kernel @ contrasts for kernel in kernels]
    raw = np.asarray([[np.sum(a*b) for b in kernels] for a in kernels])
    projected = np.asarray([[np.sum(a*b) for b in restricted] for a in restricted])
    raw_error = np.max(abs(raw-np.asarray(audit['raw_gram'])), initial=0)
    projected_error = np.max(abs(projected-np.asarray(audit['projected_gram'])), initial=0)
    np.testing.assert_allclose(audit['raw_gram'], raw, rtol=2e-12, atol=2e-10)
    np.testing.assert_allclose(audit['projected_gram'], projected, rtol=3e-11, atol=2e-8)
    assert np.all(abs(projected-np.asarray(audit['projected_gram'])) <=
                  np.asarray(audit['projected_roundoff_envelope']) + 1e-10)
    return audit, float(raw_error), float(projected_error), float(
        np.linalg.norm(projected-np.asarray(audit['projected_gram'])) /
        max(np.linalg.norm(projected), np.finfo(float).tiny))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(20261003); n = 18; errors = []
    def check(*values):
        audit, raw_error, projected_error, relative_error = dense_check(*values)
        errors.append((raw_error, projected_error, relative_error))
        return audit
    for entity_scale in [.01, 1., 100.]:
        for rank in [0, 1, 4]:
            labels = np.repeat(np.arange(3), 6)
            z = np.zeros((n, 12)); w = np.zeros((n, 9))
            for block in range(3):
                z[block*6:(block+1)*6, block*4:(block+1)*4] = rng.normal(size=(6, 4))*entity_scale
                w[block*6:(block+1)*6, block*3:(block+1)*3] = rng.normal(size=(6, 3))
            # Duplicate COO cells cancel and coalesce without changing the kernel.
            rr, cc = np.nonzero(z)
            duplicate = sparse.coo_matrix((np.r_[z[rr,cc], np.ones(len(rr)), -np.ones(len(rr))],
                                          (np.tile(rr,3), np.tile(cc,3))), shape=z.shape)
            operators = {'signed': duplicate, 'reuse': sparse.csr_matrix(w),
                         'zero': sparse.csr_matrix((n, 5))}
            factor = rng.normal(size=(n, rank)); diagonal = rng.uniform(.3, 2., n)
            x = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
            for permutation in [np.arange(n), rng.permutation(n)]:
                selected = {name: value.tocsr()[permutation] for name,value in operators.items()}
                audit = check(labels[permutation], selected,
                    factor[permutation], diagonal[permutation], x[permutation])
                assert audit['exactly_zero_incidence_names'] == ['zero']
                assert not audit['covariance_model_selected']
    labels=np.zeros(n,dtype=int); x=np.ones((n,1)); f=rng.normal(size=(n,2)); d=np.ones(n)
    target={'target':sparse.eye(n,format='csr')}
    uniform=check(labels,target,f,d,x)
    assert uniform['raw_diagnostics']['rank']==2
    assert uniform['reml_diagnostics']['rank']==2
    weighted=check(labels,target,f,np.linspace(.5,2.,n),x)
    assert weighted['raw_diagnostics']['rank']==3
    assert weighted['reml_diagnostics']['rank']==3
    # A nonzero raw family kernel is annihilated by the fixed intercept.
    intercept=check(labels,{'family':sparse.csr_matrix(np.ones((n,1)))},f,d,x)
    assert intercept['raw_diagnostics']['rank']==3
    assert intercept['reml_diagnostics']['unresolved_kernel_names']==['family']
    assert intercept['reml_diagnostics']['disposition']=='unresolved_kernel_norm_requires_review'
    # Independent raw bases become proportional after removing fixed effects.
    mixed=sparse.hstack([sparse.eye(n),sparse.csr_matrix(np.ones((n,1)))],format='csr')
    projected=check(labels,{'identity_plus_intercept':mixed},f,d,x)
    assert projected['raw_diagnostics']['rank']==3
    assert projected['reml_diagnostics']['rank']==2
    assert projected['reml_diagnostics']['disposition']=='dependent_covariance_bases_require_review'
    # Equivalent column spaces, despite column scaling, preserve both Grams.
    design=np.column_stack([np.ones(n),rng.normal(size=n)])
    base=check(labels,target,f,d,design)
    scaled=check(labels,target,f,d,design*np.array([1e200,1e-200]))
    np.testing.assert_allclose(base['projected_gram'],scaled['projected_gram'],rtol=2e-12,atol=1e-10)
    dense_shapes=[]; original_toarray=sparse.csr_matrix.toarray
    def tracked_toarray(matrix, *args, **kwargs):
        dense_shapes.append(matrix.shape)
        return original_toarray(matrix,*args,**kwargs)
    with patch.object(sparse.csr_matrix,'toarray',tracked_toarray):
        audit_covariance_bases(np.repeat(np.arange(3),6),target,f,d,x)
    assert dense_shapes and all(shape[0]<=6 and shape[1]<=6 for shape in dense_shapes)
    invalid=[]
    def rejected(name, **replacement):
        values=dict(block_labels=labels,incidence=target,species_factor=f,residual_diagonal=d,design=x)
        values.update(replacement)
        try: audit_covariance_bases(**values)
        except (ValueError,ArithmeticError): invalid.append(name)
        else: raise AssertionError('Invalid input accepted: '+name)
    rejected('empty_rows',block_labels=[])
    rejected('nonpositive_diagonal',residual_diagonal=np.zeros(n))
    rejected('nonfinite_diagonal',residual_diagonal=np.full(n,np.inf))
    rejected('factor_shape',species_factor=np.ones((n-1,2)))
    rejected('nonfinite_factor',species_factor=np.full((n,1),np.nan))
    rejected('zero_design_column',design=np.zeros((n,1)))
    rejected('rank_deficient_design',design=np.ones((n,2)))
    rejected('nonfinite_design',design=np.full((n,1),np.nan))
    rejected('no_residual_dimension',design=np.eye(n))
    rejected('reserved_entity_name',incidence={'species':sparse.eye(n)})
    rejected('wrong_entity_rows',incidence={'a':sparse.eye(n-1)})
    rejected('nonfinite_incidence',incidence={'a':sparse.csr_matrix(np.full((n,1),np.inf))})
    rejected('cross_component_entity',block_labels=np.arange(n),incidence={'a':sparse.csr_matrix(np.ones((n,1)))})
    receipt=dict(status='passed_covariance_basis_error_contrast_checks',
        checked_utc=datetime.now(timezone.utc).isoformat(),dense_cases=len(errors),
        maximum_raw_gram_absolute_error=max(e[0] for e in errors),
        maximum_projected_gram_absolute_error=max(e[1] for e in errors),
        maximum_projected_gram_relative_frobenius_error=max(e[2] for e in errors),
        largest_component_dense_shape_checked=list(max(dense_shapes)),
        rejected_invalid_inputs=invalid,
        scientific_regressions=['target_identity_aliases_uniform_residual',
          'nonuniform_residual_requires_new_qualification','family_intercept_annihilated_by_fixed_effect',
          'raw_independent_kernels_alias_in_reml','extreme_design_scaling_preserves_column_space'],
        source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(__file__),Path(__file__).with_name('covariance_basis_audit.py')]},
        scope='Software algebra checks against independent dense null-space error contrasts. '
              'No production cohort acceptance, covariance selection, variance fitting or pilot.')
    with args.output.open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(receipt))


if __name__=='__main__':main()
