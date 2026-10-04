"""Detect in-memory numeric mutation without changing likelihood arithmetic.

The complete cached source is checked at cohort boundaries. Each native fit
or independent numeric replay also checks its supplied rows, design, response,
operators and residual diagonal. Changed arrays and both digest inventories
are retained before a failed cohort can receive a success checkpoint.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from full_weighted_covariance_qualification_parallel_v1 import arrays, array_bytes
from full_weighted_fit_exports import atomic


def references(value, prefix):
    seen = set()
    result = []
    for name, array in arrays(value, prefix):
        if id(array) not in seen:
            seen.add(id(array))
            result.append((name, array))
    assert result, 'Numeric guard has no arrays'
    return result


def inventory(refs):
    return {name: dict(sha256=array_bytes(value), strides=list(value.strides))
            for name, value in refs}


class ArrayGuard:
    def __init__(self, value, prefix):
        self.refs = references(value, prefix)
        self.before = inventory(self.refs)

    def check(self, root, identity, phase):
        after = inventory(self.refs)
        if after == self.before:
            return len(self.refs)
        context = dict(identity=identity, phase=phase)
        tag = hashlib.sha256(json.dumps(context, sort_keys=True).encode()).hexdigest()
        destination = Path(root) / 'memory_failures' / tag
        destination.mkdir(parents=True, exist_ok=False)
        changed = [name for name in self.before if self.before[name] != after[name]]
        files = {}
        for name, value in self.refs:
            if name in changed:
                filename = hashlib.sha256(name.encode()).hexdigest() + '.npy'
                with (destination / filename).open('xb') as handle:
                    np.save(handle, value, allow_pickle=False)
                files[name] = filename
        atomic(destination / 'failure.json', dict(**context, before=self.before,
            after=after, changed_arrays=changed, captured_changed_arrays=files,
            scientific_eligibility=False,
            scope='Digests describe pre-call numeric values and layouts. NPY files '
                  'retain changed values; closed source files retain original inputs. '
                  'This detects mutation and does not repair its native cause.'))
        raise AssertionError('Weighted fit numeric input changed: ' + phase)


def guarded_call(function, source, plan, rows, identity, matrix, response,
                 operators, audit, route, diagonal, *, root, exported=None,
                 reader=False):
    guard = ArrayGuard((rows, matrix, response, operators, diagonal), 'native-input')
    try:
        if reader:
            return function(source, plan, rows, identity, matrix, response,
                            exported, operators, audit, route, diagonal)
        return function(source, plan, rows, identity, matrix, response,
                        operators, audit, route, diagonal)
    finally:
        guard.check(root, identity, 'independent-native-after' if reader else 'producer-native-after')
