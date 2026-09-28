#!/usr/bin/env python3
"""Read every saved alignment against the same-process tree; no convergence claim."""
import csv
from io import StringIO
import json
import math
from pathlib import Path

from Bio import Phylo, SeqIO
from ancestral_chain_attempt import sha
from audit_baliphy_sample_mapping import descendant_index


def blocks(path):
    iteration, lines = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('iterations = '):
                if iteration is not None:
                    yield iteration, ''.join(lines)
                iteration, lines = int(line.split('=', 1)[1]), []
            else:
                assert iteration is not None or not line.strip()
                lines.append(line)
    if iteration is not None:
        yield iteration, ''.join(lines)


def readback(chain, receipt_path, iterations, mapping_path):
    receipt_path = Path(receipt_path)
    receipt = json.loads(receipt_path.read_text())
    assert receipt['exit_code'] == 0
    for name, digest in receipt['artifacts'].items():
        assert sha(receipt_path.parent / name) == digest
    for field in ['tree', 'alignment']:
        assert sha(chain[field]) == chain[field + '_sha256']
    directories = list(receipt_path.parent.glob('independent-chain-*'))
    assert len(directories) == 1
    directory = directories[0]
    source = Phylo.read(chain['tree'], 'newick')
    runtime = Phylo.read(directory / 'runtime-tree.nwk', 'newick')
    si, ri = descendant_index(source), descendant_index(runtime)
    assert set(si) == set(ri)
    for key, node in si.items():
        if node is not source.root:
            assert abs(node.branch_length - ri[key].branch_length) < 1e-10
    labels = {node.name for node in runtime.find_clades()}
    assert len(labels) == len(ri)
    records = list(SeqIO.parse(chain['alignment'], 'fasta'))
    observed = {r.id: str(r.seq).replace('-', '').upper() for r in records}
    assert len(records) == len(observed) == chain['proteins']
    assert set(observed) == {node.name for node in runtime.get_terminals()}
    dataset = 'whole' if '-whole-' in chain['original_configuration_ids'][0] else 'domain'
    mappings = [r for r in csv.DictReader(Path(mapping_path).open(), delimiter='\t')
                if r['guide'] == 'profile' and r['family'] == chain['family'] and r['dataset'] == dataset]
    assert len(mappings) == 4
    with (directory / 'C1.log').open() as handle:
        logs = list(csv.DictReader(handle, delimiter='\t'))
    assert [int(r['iter']) for r in logs] == list(range(iterations + 1))
    for row in logs:
        scores = [float(row[k]) for k in ['prior', 'likelihood', 'posterior']]
        assert all(math.isfinite(v) for v in scores)
        assert abs(scores[0] + scores[1] - scores[2]) < 1e-7
    sample_rows, seen = [], []
    for iteration, text in blocks(directory / 'C1.P1.fastas'):
        seen.append(iteration)
        records = list(SeqIO.parse(StringIO(text), 'fasta'))
        sequences = {r.id: str(r.seq).upper() for r in records}
        assert len(records) == len(sequences) == len(labels) and set(sequences) == labels
        assert len({len(s) for s in sequences.values()}) == 1
        assert all(set(s) <= set('ARNDCQEGHILKMFPSTWYVX-') for s in sequences.values())
        for name, expected in observed.items():
            actual = sequences[name].replace('-', '')
            assert len(actual) == len(expected)
            assert all(a == b or a == 'X' for a, b in zip(expected, actual))
        for mapping in mappings:
            key = tuple(sorted(json.loads(mapping['retained_set_json'])))
            node = ri[key]
            sequence = sequences[node.name].replace('-', '')
            sample_rows.append(dict(iteration=iteration, level=mapping['level'],
                source_node=mapping['source_node'], runtime_node=node.name,
                ungapped_length=len(sequence)))
    assert seen == list(range(0, iterations + 1, 10))
    return dict(status='all_saved_alignments_and_candidate_nodes_checked',
        attempt_receipt=str(receipt_path), attempt_receipt_sha256=sha(receipt_path),
        scalar_log=str(directory / 'C1.log'), scalar_log_sha256=sha(directory / 'C1.log'),
        candidate_samples=sample_rows, iterations=iterations,
        mapping_sha256=sha(mapping_path),
        scope='Integrity and node identities only. No stationary posterior or convergence qualification.')
