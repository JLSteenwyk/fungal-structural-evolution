"""Separate native FASTA/Newick decoding and extant-residue projection.

No Bio.SeqIO, original block parser, anchor builder or projection imports.
Source-to-runtime candidate mapping is inherited from the closed audit.
"""
import math
import re

import numpy as np

ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX-'
_ITERATION = re.compile(r'iterations\s*=\s*(0|[1-9][0-9]*)\Z')


def fasta_records(lines):
    records = {}; label = None; pieces = []
    def finish():
        if label is not None:
            sequence = ''.join(pieces).upper()
            if not sequence or any(c not in ALPHABET for c in sequence):
                raise ValueError('Empty or invalid sequence: ' + label)
            records[label] = sequence
    for raw in lines:
        line = raw.strip()
        if not line: continue
        if line.startswith('>'):
            finish(); header = line[1:].split()
            if not header or header[0] in records:
                raise ValueError('Empty or duplicate FASTA identifier')
            label = header[0]; pieces = []
        else:
            if label is None: raise ValueError('Sequence before FASTA header')
            pieces.append(''.join(line.split()))
    finish()
    if not records: raise ValueError('Empty FASTA')
    return records


def native_alignments(lines):
    iteration = None; pending = []; last = -1
    for raw in lines:
        match = _ITERATION.fullmatch(raw.strip())
        if match:
            current = int(match.group(1))
            if current <= last: raise ValueError('Repeated or unordered native iteration')
            if iteration is not None: yield iteration, fasta_records(pending)
            iteration = current; last = current; pending = []
        elif iteration is None:
            if raw.strip(): raise ValueError('Unexpected native preamble')
        else: pending.append(raw)
    if iteration is None: raise ValueError('Missing native iteration marker')
    yield iteration, fasta_records(pending)


def newick_tokens(text):
    i = 0
    while i < len(text):
        if text[i].isspace(): i += 1; continue
        if text[i] == '[':
            depth = 1; i += 1
            while i < len(text) and depth:
                depth += (text[i] == '[') - (text[i] == ']'); i += 1
            if depth: raise ValueError('Unclosed Newick comment')
        elif text[i] in '(),:;':
            yield text[i]; i += 1
        elif text[i] == "'":
            pieces = []; i += 1
            while True:
                if i == len(text): raise ValueError('Unclosed quoted Newick label')
                if text[i] == "'":
                    if i + 1 < len(text) and text[i + 1] == "'": pieces.append("'"); i += 2
                    else: i += 1; break
                else: pieces.append(text[i]); i += 1
            yield ('label', ''.join(pieces))
        else:
            start = i
            while i < len(text) and not text[i].isspace() and text[i] not in '(),:;[]': i += 1
            if start == i: raise ValueError('Invalid Newick token')
            yield ('label', text[start:i])


def tree_labels(text):
    """Iterative grammar retains every named node, including multifurcations."""
    tokens = list(newick_tokens(text)); index = 0; stack = []; nodes = []; root = None
    expect_child = True
    while index < len(tokens):
        token = tokens[index]
        if expect_child:
            if token == '(':
                stack.append([]); index += 1; continue
            if not isinstance(token, tuple): raise ValueError('Expected Newick subtree')
            children = []; label = token[1]; index += 1
        else:
            if token == ',':
                if not stack: raise ValueError('Comma outside Newick clade')
                expect_child = True; index += 1; continue
            if token == ')':
                if not stack or not stack[-1]: raise ValueError('Empty Newick clade')
                children = stack.pop(); index += 1
                if index >= len(tokens) or not isinstance(tokens[index], tuple):
                    raise ValueError('Native runtime node must have a label')
                label = tokens[index][1]; index += 1
            elif token == ';':
                if stack or root is None or index != len(tokens) - 1:
                    raise ValueError('Incomplete or extra Newick tree')
                break
            else: raise ValueError('Unexpected Newick suffix')
        if not label: raise ValueError('Empty runtime label')
        if index < len(tokens) and tokens[index] == ':':
            index += 1
            if index >= len(tokens) or not isinstance(tokens[index], tuple): raise ValueError('Missing branch length')
            length = float(tokens[index][1]); index += 1
            if not math.isfinite(length): raise ValueError('Invalid native branch length')
        node = dict(label=label, children=children); nodes.append(node)
        if stack: stack[-1].append(node)
        else:
            if root is not None: raise ValueError('Multiple Newick roots')
            root = node
        expect_child = False
    else: raise ValueError('Newick terminator required')
    labels = [n['label'] for n in nodes]
    if len(labels) != len(set(labels)): raise ValueError('Duplicate runtime node label')
    return set(labels), {n['label'] for n in nodes if not n['children']}


def project(sequences, observed, candidates):
    if not observed or not candidates: raise ValueError('Tips and candidates required')
    if len(set(candidates.values())) != len(candidates): raise ValueError('Repeated candidate label')
    if set(observed) & set(candidates.values()): raise ValueError('Candidate overlaps tip')
    if not (set(observed) | set(candidates.values())) <= set(sequences): raise ValueError('Missing tip/candidate')
    widths = {len(s) for s in sequences.values()}
    if len(widths) != 1 or not all(s and set(s) <= set(ALPHABET) for s in sequences.values()):
        raise ValueError('Invalid saved alignment')
    width = widths.pop(); occupied = np.zeros(width, dtype=bool); positions = []
    for tip in sorted(observed):
        expected = observed[tip]
        if not expected or '-' in expected or not set(expected) <= set(ALPHABET):
            raise ValueError('Nonempty ungapped observed sequence required')
        actual = sequences[tip].replace('-', '')
        if len(actual) != len(expected) or any(a != b and a != 'X' for a,b in zip(expected, actual)):
            raise ValueError('Changed observed residue identity')
        row = np.frombuffer(sequences[tip].encode('ascii'), dtype=np.uint8)
        selected = np.flatnonzero(row != ord('-')); positions.append(selected); occupied[selected] = True
    positions = np.concatenate(positions)
    lookup = np.full(256, 255, dtype=np.uint8)
    lookup[np.frombuffer(ALPHABET.encode(), dtype=np.uint8)] = np.arange(len(ALPHABET), dtype=np.uint8)
    rows = np.stack([np.frombuffer(sequences[candidates[node]].encode(), dtype=np.uint8) for node in sorted(candidates)])
    states = lookup[rows[:, positions]]
    unanchored = np.count_nonzero((rows != ord('-')) & ~occupied, axis=1).astype(np.int32)
    assert np.all(states < len(ALPHABET))
    return states, unanchored
