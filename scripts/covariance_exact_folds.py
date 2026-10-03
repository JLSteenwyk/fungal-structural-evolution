"""Exact dyadic kernel certificates and nonnegative uniform variance maps.

Named identities authorize folds; numerical ranks never authorize deletion.
The global operator files remain unchanged.
"""
import numpy as np
from scipy import sparse

ORIGINAL_NAMES = ['residual', 'target_node', 'background_node', 'model_pair',
                  'gene', 'model', 'family', 'family_intercept', 'species']
RELATIONS = ['pair_equals_target_plus_background', 'gene_equals_half_pair', 'model_equals_gene']


def integer_operator(value, scale):
    value = sparse.csr_matrix(value, copy=True)
    assert value.has_canonical_format and np.isfinite(value.data).all()
    scaled = value.data * scale
    assert np.array_equal(scaled, np.rint(scaled)) and np.max(np.abs(scaled), initial=0) <= 2
    assert np.max(np.diff(value.indptr), initial=0) <= 4
    return sparse.csr_matrix((scaled.astype('int64'), value.indices.copy(), value.indptr.copy()), shape=value.shape)


def difference_record(matrix, rows):
    matrix = matrix.tocsr(); matrix.eliminate_zeros(); matrix.sort_indices()
    if matrix.nnz == 0:
        return dict(exact=True, nonzero_entries=0, witness=None)
    row = int(np.flatnonzero(np.diff(matrix.indptr))[0]); start = matrix.indptr[row]
    column = int(matrix.indices[start])
    return dict(exact=False, nonzero_entries=int(matrix.nnz), witness=dict(
        row=row, column=column, original_case_row=int(rows[row]), original_case_column=int(rows[column]),
        scaled_integer_difference=int(matrix.data[start])))


def certificate(operators, family_intercept, mode, rows):
    assert mode in ['signed', 'unsigned']
    n = len(rows); assert n > 0
    integer = {name: integer_operator(operators[name], 2 if name in ['gene', 'model'] else 1)
               for name in ['target_node', 'background_node', 'model_pair', 'gene', 'model']}
    assert all(value.shape[0] == n for value in integer.values())
    target = integer['target_node']
    assert np.all(np.diff(target.indptr) == 1) and np.all(target.data == 1)
    assert len(np.unique(target.indices)) == n
    family = integer_operator(operators['family'], 1)
    intercept = integer_operator(family_intercept, 1)
    assert family.shape == intercept.shape and np.all(intercept.data == 1)
    assert np.all(np.diff(intercept.indptr) == 1)
    assert family.nnz == 0 if mode == 'signed' else (family - 2 * intercept).nnz == 0
    kernels = {name: value @ value.T for name, value in integer.items()}
    checks = {
        RELATIONS[0]: difference_record(kernels['model_pair'] - kernels['target_node'] - kernels['background_node'], rows),
        RELATIONS[1]: difference_record(kernels['gene'] - 2 * kernels['model_pair'], rows),
        RELATIONS[2]: difference_record(kernels['model'] - kernels['gene'], rows)}
    return dict(records=n, loading_mode=mode, target_identity_exact=True,
        family_fold_exact=True, relations=checks,
        integer_scales={'target_node': 1, 'background_node': 1, 'model_pair': 1, 'gene': 2, 'model': 2},
        interpretation='Exact integer Gram equality; no numerical tolerance, rank-based deletion or estimated covariance.')


def variance_map(checks, mode):
    assert mode in ['signed', 'unsigned']
    active = ['residual', 'background_node', 'model_pair', 'gene', 'model', 'family_intercept', 'species']
    coefficients = {name: np.eye(9, dtype='float64')[ORIGINAL_NAMES.index(name)].copy() for name in active}
    coefficients['residual'][ORIGINAL_NAMES.index('target_node')] = 1
    coefficients['family_intercept'][ORIGINAL_NAMES.index('family')] = 0 if mode == 'signed' else 4
    folds = []
    # Process dependencies in reverse order so later composite columns keep
    # every earlier original contribution. Every coefficient is exactly dyadic.
    for relation, removed, targets in [
        (RELATIONS[2], 'model', [('gene', 1)]),
        (RELATIONS[1], 'gene', [('model_pair', .5)]),
        (RELATIONS[0], 'model_pair', [('residual', 1), ('background_node', 1)])]:
        assert type(checks[relation]['exact']) is bool
        if checks[relation]['exact']:
            value = coefficients.pop(removed); active.remove(removed)
            for name, amount in targets: coefficients[name] += amount * value
            folds.append(dict(identity=relation, removed=removed, targets=dict(targets)))
    transform = np.asarray([coefficients[name] for name in active])
    inverse = np.zeros((9, len(active)))
    for index, name in enumerate(active): inverse[ORIGINAL_NAMES.index(name), index] = 1
    assert np.all(transform >= 0) and np.all(inverse >= 0)
    assert np.array_equal(transform @ inverse, np.eye(len(active)))
    return dict(original_names=ORIGINAL_NAMES, retained_names=active, forward=transform.tolist(),
        nonnegative_right_inverse=inverse.tolist(), folds=folds,
        covariance_cone_preserved_given_exact_certificates=True,
        component_variance_attribution_accepted=False, nonuniform_weighting_accepted=False)


def independent_gram(value, scale):
    """Exact independent CSC column outer sums, with vectorized singletons."""
    matrix = sparse.csc_matrix(value, copy=True)
    assert matrix.has_canonical_format and np.isfinite(matrix.data).all()
    weights = matrix.data * scale
    assert np.array_equal(weights, np.rint(weights)) and np.max(abs(weights), initial=0) <= 2
    weights = weights.astype('int64'); sizes = np.diff(matrix.indptr)
    singleton = matrix.indptr[:-1][sizes == 1]
    diagonal = np.zeros(matrix.shape[0], dtype='int64')
    np.add.at(diagonal, matrix.indices[singleton], weights[singleton] ** 2)
    row_arrays = [np.arange(matrix.shape[0])]; column_arrays = [np.arange(matrix.shape[0])]; data_arrays = [diagonal]
    for column in np.flatnonzero(sizes > 1):
        a, b = matrix.indptr[column:column + 2]
        rows = matrix.indices[a:b]; values = weights[a:b]
        row_arrays.append(np.repeat(rows, len(rows)))
        column_arrays.append(np.tile(rows, len(rows)))
        data_arrays.append((values[:, None] * values[None, :]).ravel())
    result = sparse.coo_matrix((np.concatenate(data_arrays),
        (np.concatenate(row_arrays), np.concatenate(column_arrays))), shape=(matrix.shape[0], matrix.shape[0]), dtype='int64').tocsr()
    result.eliminate_zeros(); result.sort_indices()
    return result


def independent_certificate(operators, family_intercept, mode, rows):
    n = len(rows)
    target = sparse.csr_matrix(operators['target_node'])
    assert target.shape[0] == n and target.nnz == n
    assert np.all(np.diff(target.indptr) == 1) and set(target.data) == {1}
    assert len(set(target.indices.tolist())) == n
    family = sparse.csr_matrix(operators['family']); intercept = sparse.csr_matrix(family_intercept)
    assert intercept.shape[0] == n and np.all(np.diff(intercept.indptr) == 1) and set(intercept.data) == {1}
    if mode == 'signed': assert family.nnz == 0
    else:
        assert mode == 'unsigned' and family.shape == intercept.shape
        assert np.array_equal(family.indptr, intercept.indptr) and np.array_equal(family.indices, intercept.indices)
        assert np.array_equal(family.data, 2 * intercept.data)
    background = independent_gram(operators['background_node'], 1)
    pair = independent_gram(operators['model_pair'], 1)
    gene = independent_gram(operators['gene'], 2)
    model = independent_gram(operators['model'], 2)
    differences = [pair - sparse.eye(n, dtype='int64', format='csr') - background,
                   gene - 2 * pair, model - gene]
    checks = {}
    for name, difference in zip(RELATIONS, differences):
        difference = difference.tocoo(); keep = difference.data != 0
        r, c, v = difference.row[keep], difference.col[keep], difference.data[keep]
        if len(v):
            index = np.lexsort((c, r))[0]
            witness = dict(row=int(r[index]), column=int(c[index]), original_case_row=int(rows[r[index]]),
                original_case_column=int(rows[c[index]]), scaled_integer_difference=int(v[index]))
        else: witness = None
        checks[name] = dict(exact=len(v) == 0, nonzero_entries=len(v), witness=witness)
    return dict(records=n, loading_mode=mode, target_identity_exact=True, family_fold_exact=True,
        relations=checks, integer_scales={'target_node': 1, 'background_node': 1, 'model_pair': 1, 'gene': 2, 'model': 2},
        interpretation='Exact integer Gram equality; no numerical tolerance, rank-based deletion or estimated covariance.')


def independent_variance_map(checks, mode):
    assert mode in ['signed', 'unsigned']
    assert all(type(checks[name]['exact']) is bool for name in RELATIONS)
    pair = {'residual': 1, 'background_node': 1} if checks[RELATIONS[0]]['exact'] else {'model_pair': 1}
    gene = {name: .5 * value for name, value in pair.items()} if checks[RELATIONS[1]]['exact'] else {'gene': 1}
    model = dict(gene) if checks[RELATIONS[2]]['exact'] else {'model': 1}
    expressions = [{'residual': 1}, {'residual': 1}, {'background_node': 1}, pair, gene, model,
        {} if mode == 'signed' else {'family_intercept': 4}, {'family_intercept': 1}, {'species': 1}]
    order = ['residual', 'background_node', 'model_pair', 'gene', 'model', 'family_intercept', 'species']
    names = [name for name in order if any(name in expression for expression in expressions)]
    transform = [[expression.get(name, 0) for expression in expressions] for name in names]
    inverse = [[int(original == name) for name in names] for original in ORIGINAL_NAMES]
    assert np.array_equal(np.asarray(transform) @ np.asarray(inverse), np.eye(len(names)))
    folds = []
    for identity, removed, targets in [(RELATIONS[2], 'model', {'gene': 1}),
        (RELATIONS[1], 'gene', {'model_pair': .5}),
        (RELATIONS[0], 'model_pair', {'residual': 1, 'background_node': 1})]:
        if checks[identity]['exact']: folds.append(dict(identity=identity, removed=removed, targets=targets))
    return dict(original_names=ORIGINAL_NAMES, retained_names=names, forward=transform,
        nonnegative_right_inverse=inverse, folds=folds,
        covariance_cone_preserved_given_exact_certificates=True,
        component_variance_attribution_accepted=False, nonuniform_weighting_accepted=False)
