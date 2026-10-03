"""SQL partition counts and Fraction division independently reconstruct controls."""
from fractions import Fraction
import sqlite3

import numpy as np

from inverse_reuse_weight_controls import POLICIES


def database(source):
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE cases(case_row INTEGER PRIMARY KEY, bg TEXT, pair TEXT, component TEXT)')
    keys = source['keys']
    db.executemany('INSERT INTO cases VALUES(?,?,?,?)',
        ((i, keys['background_node'][i].decode(), keys['background_pair'][i].decode(),
          keys['family_component'][i].decode()) for i in range(len(source['ids']))))
    db.execute('CREATE TEMP TABLE membership(ordinal INTEGER PRIMARY KEY, case_row INTEGER UNIQUE)')
    return db


def reconstruct(db, selected):
    n = len(selected); assert n > 0
    db.execute('DELETE FROM membership')
    db.executemany('INSERT INTO membership VALUES(?,?)', ((i, int(v)) for i, v in enumerate(selected)))
    result = db.execute('''SELECT ordinal, COUNT(*) OVER(PARTITION BY bg),
        COUNT(*) OVER(PARTITION BY pair), COUNT(*) OVER(PARTITION BY component)
        FROM membership JOIN cases USING(case_row) ORDER BY ordinal''').fetchall()
    assert [r[0] for r in result] == list(range(n))
    counts = np.asarray([r[1:] for r in result], dtype=np.int64).T
    groups = list(db.execute('''SELECT COUNT(DISTINCT bg), COUNT(DISTINCT pair),
        COUNT(DISTINCT component) FROM membership JOIN cases USING(case_row)''').fetchone())
    weights = np.empty((4, n)); diagonals = np.empty((4, n)); weights[0] = 1; diagonals[0] = 1
    for j, g in enumerate(groups):
        # Compute each distinct count via exact rationals; no producer NumPy grouping/division.
        w = {int(c): float(Fraction(n, g * int(c))) for c in set(counts[j].tolist())}
        d = {int(c): float(Fraction(g * int(c), n)) for c in set(counts[j].tolist())}
        weights[j + 1] = [w[int(c)] for c in counts[j]]
        diagonals[j + 1] = [d[int(c)] for c in counts[j]]
        assert sum(Fraction(1, int(c)) for c in counts[j]) == g
    return counts, weights, diagonals, groups


def check_policy_records(saved, n, counts, weights, diagonals, groups):
    assert len(saved) == 4
    for j, r in enumerate(saved):
        reuse = [1] * n if j == 0 else counts[j - 1].tolist()
        g = n if j == 0 else groups[j - 1]
        one = all(g * c == n for c in reuse)
        expected = dict(policy=POLICIES[j], represented_groups=g, records=n,
            maximum_reuse=max(reuse), minimum_reuse=min(reuse),
            weight_min=float(min(weights[j])), weight_max=float(max(weights[j])),
            diagonal_min=float(min(diagonals[j])), diagonal_max=float(max(diagonals[j])),
            weight_sum_float64=float(np.sum(weights[j])), exact_weight_sum=n,
            exact_group_total_numerator=n, exact_group_total_denominator=g,
            diagonal_is_exact_uniform_one=one,
            future_basis_disposition=('requires_uniform_named_fold_and_fresh_numerical_qualification'
                if one else 'requires_positive_diagonal_cone_and_fresh_numerical_qualification'))
        assert r == expected
        assert abs(float(np.sum(weights[j])) - n) <= 16 * np.finfo(float).eps * n
        assert np.all(weights[j] > 0) and np.all(diagonals[j] > 0)
        assert np.max(np.abs(weights[j] * diagonals[j] - 1)) <= 4 * np.finfo(float).eps
