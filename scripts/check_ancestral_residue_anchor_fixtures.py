#!/usr/bin/env python3
"""Validate anchor projection on every frozen effective-input short-run fixture."""
import argparse
import csv
from io import StringIO
import json
from pathlib import Path
from Bio import Phylo, SeqIO
from ancestral_chain_attempt import sha, write_json
from ancestral_residue_anchors import anchored_columns
from audit_baliphy_sample_mapping import descendant_index
from readback_independent_baliphy_chain import blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    inputs = Path('results/ancestral/baliphy-independent-chain-inputs-20260927-v1/chain_inputs.json')
    proof = json.loads(Path('metadata/baliphy_sample_mapping_final_readback_completed_20260927.json').read_text())
    ap = Path(proof['audit_receipt'])
    assert sha(ap) == proof['audit_receipt_sha256']
    audited = json.loads(ap.read_text())
    mapping = Path('results/ancestral/case-local-trees-20260927-v1/ancestral_node_mapping.tsv')
    with mapping.open() as handle:
        mappings = list(csv.DictReader(handle, delimiter='\t'))
    selected = {}
    for row in json.loads(inputs.read_text()):
        selected.setdefault(row['effective_input_group'], row)
    assert len(selected) == 135
    counts = dict(inputs=0, alignments=0, candidate_reconstructions=0,
                  anchored_columns=0, unanchored_columns=0, extant_residue_observations=0)
    pins = {str(inputs): sha(inputs), str(ap): sha(ap), str(mapping): sha(mapping)}
    for row in selected.values():
        root = Path('results/ancestral/baliphy-sample-mapping-20260927-v1') / row['original_configuration_ids'][0]
        rp = root / 'receipt.json'
        assert sha(rp) == audited['pins'][str(rp)]
        pins[str(rp)] = sha(rp)
        receipt = json.loads(rp.read_text())
        for name, digest in receipt['artifacts'].items():
            assert sha(root / name) == digest
        assert sha(row['alignment']) == row['alignment_sha256']
        runtime = descendant_index(Phylo.read(root / 'runtime-tree.nwk', 'newick'))
        dataset = 'whole' if '-whole-' in row['original_configuration_ids'][0] else 'domain'
        chosen = [m for m in mappings if m['guide'] == 'profile'
                  and m['family'] == row['family'] and m['dataset'] == dataset]
        assert len(chosen) == 4
        candidates = {m['source_node']: runtime[tuple(sorted(json.loads(m['retained_set_json'])))].name
                      for m in chosen}
        assert len(candidates) == 4
        observed = {r.id: str(r.seq).replace('-', '').upper()
                    for r in SeqIO.parse(row['alignment'], 'fasta')}
        samples = list(root.glob('mapping-check-*/C1.P1.fastas'))
        assert len(samples) == 1
        iterations = []
        for iteration, text in blocks(samples[0]):
            iterations.append(iteration)
            records = list(SeqIO.parse(StringIO(text.lstrip()), 'fasta'))
            sequences = {r.id: str(r.seq).upper() for r in records}
            assert len(sequences) == len(records)
            columns = anchored_columns(sequences, observed, candidates)
            # Independent traversal verifies that every input residue appears
            # exactly once, in order; column shifts cannot lose/duplicate it.
            coordinates = {tip: [] for tip in observed}
            for column in columns:
                for tip, position in column['anchors']:
                    coordinates[tip].append(position)
                counts['unanchored_columns' if column['unanchored'] else 'anchored_columns'] += 1
            for tip, sequence in observed.items():
                assert coordinates[tip] == list(range(1, len(sequence) + 1))
                counts['extant_residue_observations'] += len(sequence)
            for node, label in candidates.items():
                actual = ''.join(c['states'][node] for c in columns).replace('-', '')
                assert actual == sequences[label].replace('-', '')
                counts['candidate_reconstructions'] += 1
            counts['alignments'] += 1
        assert iterations == [0, 10, 20]
        counts['inputs'] += 1
    for script in [__file__, 'scripts/ancestral_residue_anchors.py',
                   'scripts/test_ancestral_residue_anchors.py']:
        pins[script] = sha(script)
    write_json(args.output, dict(status='passed_extant_anchor_projection_fixtures', counts=counts,
        pins=pins, scope='All 135 effective short-run inputs; residue coverage/order and four candidate '
        'sequence reconstructions checked for every saved alignment. Coordinate provenance is input '
        'ungapped sequence, not whole-protein coordinates. No production-chain mixing diagnostics, '
        'cross-sample identities for unanchored columns, or posterior qualification.'))
    print(json.dumps(counts, indent=2))


if __name__ == '__main__':
    main()
