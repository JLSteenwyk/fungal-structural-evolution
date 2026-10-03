#!/usr/bin/env python3
"""Independent decoder/projection contracts against the original implementation."""
import argparse
from io import StringIO
import json
from pathlib import Path
import tempfile

from Bio import Phylo, SeqIO
import numpy as np

from independent_native_ancestral_alignment import ALPHABET, fasta_records, native_alignments, tree_labels, project
from prepare_ancestral_state_traces import project_states
from readback_independent_baliphy_chain import blocks
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); rng = np.random.default_rng(20261005); cases = 0
    for tips in [1, 2, 5, 11]:
        for width in [1, 3, 17, 61]:
            for repeat in range(5):
                labels = ['tip-' + str(i) for i in range(tips)] + ['runtime-' + str(i) for i in range(4)]
                matrix = rng.integers(0, len(ALPHABET), (len(labels), width))
                matrix[:tips, 0] %= len(ALPHABET) - 1
                sequences = {label:''.join(ALPHABET[s] for s in row) for label,row in zip(labels,matrix)}
                observed = {tip:sequence.replace('-', '') for tip,sequence in sequences.items() if tip.startswith('tip-')}
                if repeat == 1: observed = {tip:'X' * len(seq) for tip,seq in observed.items()}
                candidates = {'source-' + str(3-i):'runtime-' + str(i) for i in range(4)}
                a,b = project_states(sequences, observed, candidates); c,d = project(sequences, observed, candidates)
                assert np.array_equal(a,c) and np.array_equal(b,d)
                text = '\n'.join('>' + label + ' description\n' + '\n'.join(sequence[j:j+7].lower() for j in range(0,len(sequence),7))
                    for label,sequence in reversed(list(sequences.items()))) + '\n'
                parsed = fasta_records(StringIO(text))
                original = {r.id:str(r.seq).upper() for r in SeqIO.parse(StringIO(text), 'fasta')}
                assert parsed == original == sequences
                c,d = project(parsed, observed, candidates); assert np.array_equal(a,c) and np.array_equal(b,d)
                cases += 1
    # Unknown expected tips are a wildcard; concrete input residues are not.
    seq = {'t':'A--', 'n':'AX-'}
    a,b = project(seq, {'t':'X'}, {'s':'n'}); assert a.tolist() == [[0]] and b.tolist() == [1]
    seq = {'t':'-A-', 'n':'C-X'}
    a,b = project(seq, {'t':'A'}, {'s':'n'}); assert a.tolist() == [[21]] and b.tolist() == [2]
    trees = ['(t1:1,t2:2)root;', "('t,1':1e-7,'t''2':0)root:2;", '(t:0)root;',
        '((a:1,b:2)i:3,c:4,d:5)root;', '(a[comment]:1,b:2)root;']
    for text in trees:
        labels, tips = tree_labels(text); original = Phylo.read(StringIO(text), 'newick')
        assert labels == {n.name for n in original.find_clades()}
        assert tips == {n.name for n in original.get_terminals()}
    # Deep native trees must not depend on Python recursion depth.
    text = 't0:1'
    for i in range(1,1200): text = '(' + text + ',t' + str(i) + ':1)i' + str(i) + ':1'
    labels,tips = tree_labels(text + ';'); assert len(labels) == 2399 and len(tips) == 1200
    with tempfile.TemporaryDirectory(prefix='independent-native-parser-') as directory:
        p = Path(directory) / 'native.fastas'
        p.write_text('\niterations = 0\n\n>t\nA-X\n>n\nAC-\n\niterations = 10\n>n\nA--\n>t\nA-X\n')
        independent = list(native_alignments(p.open()))
        original = [(i,{r.id:str(r.seq).upper() for r in SeqIO.parse(StringIO(text), 'fasta')}) for i,text in blocks(p)]
        assert independent == original and [i for i,_ in independent] == [0,10]
    failures = []
    for text in ['A\n', '>x\nA\n>x\nA\n', '>\nA\n', '>x\n', '>x\n?\n']:
        try: fasta_records(StringIO(text))
        except ValueError: failures.append('invalid-fasta')
        else: raise AssertionError('Invalid FASTA accepted')
    for text in ['>x\nA\n', 'iterations = 0\n>x\nA\niterations = 0\n>x\nA\n',
        'iterations = 10\n>x\nA\niterations = 0\n>x\nA\n', 'iterations = -1\n>x\nA\n']:
        try: list(native_alignments(StringIO(text)))
        except ValueError: failures.append('invalid-iteration-envelope')
        else: raise AssertionError('Invalid native envelope accepted')
    for text in ['(a:1,b:2);','(a:1,a:2)r;','()r;','(a:1,b:2)r', '(a:nan,b:2)r;', '(a:1,b:2)r;extra', '(a:1,b:2)r[open;']:
        try: tree_labels(text)
        except ValueError: failures.append('invalid-runtime-tree')
        else: raise AssertionError(('Invalid native tree accepted',text))
    base = {'t':'A-', 'n':'-C'}
    bad = [(base, {'t':'C'},{'s':'n'}), (base,{'t':''},{'s':'n'}), (base,{'t':'A'},{'s':'t'}),
        (base,{'t':'A'},{'s':'missing'}), ({'t':'A','n':'CC'},{'t':'A'},{'s':'n'}),
        (base,{'t':'A'},{'s':'n','other':'n'}), ({'t':'A-','n':'?- '},{'t':'A'},{'s':'n'})]
    for values,observed,candidates in bad:
        try: project(values, observed, candidates)
        except ValueError: failures.append('invalid-projection')
        else: raise AssertionError('Invalid projection accepted')
    paths = [__file__, 'scripts/independent_native_ancestral_alignment.py', 'scripts/prepare_ancestral_state_traces.py',
        'scripts/ancestral_residue_anchors.py', 'scripts/readback_independent_baliphy_chain.py']
    result = dict(status='passed_independent_native_ancestral_alignment_contracts', projection_fixtures=cases,
        original_projection_and_fasta_parsing_agree=True, unknown_wildcard_and_unanchored_fixtures=2,
        newick_fixtures=len(trees), deep_iterative_tree_nodes=2399, native_iteration_fixture_samples=2,
        malformed_inputs_rejected=len(failures), rejection_categories=failures,
        source_hashes={str(p):sha(p) for p in paths}, scientific_eligibility=False,
        scope='Software decoding/projection contracts only. No production native alignments replayed; '
              'candidate clade mapping, sampling stationarity and biological model qualification are not established.')
    with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__': main()
