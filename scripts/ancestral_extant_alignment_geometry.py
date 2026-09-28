#!/usr/bin/env python3
"""Extant alignment geometry with stable input residue positions, including X."""
import numpy as np


def project(sequences, observed):
    """Validate observed residues; encode geometry as A/gap, not biological sequence.

An input X retains its residue position regardless of sampled amino acid.
Internal-only columns are removed. Never use this encoding for sequence fits.
"""
    if not observed or not set(observed) <= set(sequences):
        raise ValueError('Missing observed tips')
    selected = {t: sequences[t].upper() for t in sorted(observed)}
    if len({len(s) for s in selected.values()}) != 1:
        raise ValueError('Unequal aligned lengths')
    for t, s in selected.items():
        original = observed[t].upper()
        actual = s.replace('-', '')
        if not original or set(original) - set('ACDEFGHIKLMNPQRSTVWYX') or set(s) - set('ACDEFGHIKLMNPQRSTVWYX-'):
            raise ValueError('Invalid extant sequence alphabet')
        if len(actual) != len(original) or any(a != b and a != 'X' for a, b in zip(original, actual)):
            raise ValueError('Observed residue identity changed')
    matrix = np.array([list(s) for s in selected.values()]) != '-'
    matrix = matrix[:, matrix.any(axis=0)]
    return {t: ''.join(np.where(row, 'A', '-')) for t, row in zip(selected, matrix)}


def signature(encoded):
    """Row (tip, ungapped position), column target tip: partner position or zero."""
    if not encoded or any(set(s) - {'A', '-'} for s in encoded.values()):
        raise ValueError('Expected geometry-only encoding')
    if len({len(s) for s in encoded.values()}) != 1:
        raise ValueError('Unequal aligned lengths')
    matrix = np.array([list(encoded[t]) for t in sorted(encoded)]) == 'A'
    positions = np.cumsum(matrix, axis=1, dtype=np.int32) * matrix
    result = np.empty((int(matrix.sum()), len(encoded)), dtype=np.int32)
    offset = 0
    for row in matrix:
        n = int(row.sum())
        result[offset:offset+n] = positions[:, row].T
        offset += n
    return result


def signature_distance(a, b):
    if a.shape != b.shape:
        raise ValueError('Mismatched signature shapes')
    # Block rows to avoid a full additional boolean matrix at large tip counts.
    return 2 * sum(int(np.count_nonzero(a[i:i+1024] != b[i:i+1024]))
                   for i in range(0, len(a), 1024))
