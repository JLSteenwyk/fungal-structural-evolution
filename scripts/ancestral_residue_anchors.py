#!/usr/bin/env python3
"""Represent sampled columns by extant residue identities, never column numbers.

Coordinates are one-based in each *input* ungapped tip sequence. They are
not necessarily whole-protein coordinates for cropped/domain inputs.
An anchor set describes a sampled homology hypothesis. Changed sets are
different hypotheses; sharing one anchor does not establish identical
homology for the remaining residues. Columns with no extant residues have
no stable cross-sample identity and remain explicitly unanchored.
"""


ALPHABET = frozenset('ARNDCQEGHILKMFPSTWYVX-')


def anchored_columns(sequences, observed, candidates):
    """Project one alignment, mapping stable source nodes to runtime labels.

Return column-order records with immutable extant anchor tuples and states
keyed by stable source-node identity. Runtime labels and column indices are
not cross-chain biological identities. Preserve unanchored columns so that
candidate ungapped sequences remain reconstructible.
"""
    if not sequences or not observed or not candidates:
        raise ValueError('Alignment, observed tips and candidate mapping required')
    if len(set(candidates.values())) != len(candidates):
        raise ValueError('Candidate runtime labels must be distinct')
    if not set(observed) <= set(sequences) or not set(candidates.values()) <= set(sequences):
        raise ValueError('Missing observed tip or candidate')
    if set(observed) & set(candidates.values()):
        raise ValueError('Candidate must be an internal node')
    lengths = {len(s) for s in sequences.values()}
    if len(lengths) != 1 or not all(set(s) <= ALPHABET for s in sequences.values()):
        raise ValueError('Invalid alignment lengths or characters')
    for tip, original in observed.items():
        if not original or not set(original) <= ALPHABET - {'-'}:
            raise ValueError('Expected ungapped nonempty tip sequence')
        actual = sequences[tip].replace('-', '')
        if len(actual) != len(original) or any(a != b and a != 'X' for a, b in zip(original, actual)):
            raise ValueError('Observed residue identity changed: ' + tip)
    tips = sorted(observed)
    positions = dict.fromkeys(tips, 0)
    result = []
    for column in range(next(iter(lengths))):
        anchors = []
        for tip in tips:
            if sequences[tip][column] != '-':
                positions[tip] += 1
                anchors.append((tip, positions[tip]))
        result.append(dict(anchors=tuple(anchors),
                           states={node: sequences[label][column]
                                   for node, label in sorted(candidates.items())},
                           unanchored=not anchors))
    return result
