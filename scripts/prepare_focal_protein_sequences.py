#!/usr/bin/env python3
"""Extract hash-verified focal and sister proteins for annotation follow-up."""
import argparse
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from inspect_focal_domain_architectures import sha


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--architectures', type=Path, required=True)
    ap.add_argument('--proteomes', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    source = json.loads(a.architectures.read_text())
    assert source['status'] == 'complete_focal_and_sister_candidate_architecture_export'
    paths = [a.architectures, Path(__file__), Path(__file__).with_name('inspect_focal_domain_architectures.py')]
    paths += [a.proteomes / (r['taxon_id'] + '.faa') for r in source['records']]
    pins = {str(p): sha(p) for p in paths}
    records = []
    for row in source['records']:
        path = a.proteomes / (row['taxon_id'] + '.faa')
        found = [r for r in SeqIO.parse(path, 'fasta') if r.id == row['protein_id']]
        assert len(found) == 1, row['tree_label']
        sequence = str(found[0].seq)
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        assert row['sequence_id'] == 'S' + digest
        records.append(dict(label=row['tree_label'], sequence=sequence, length=len(sequence),
                            sequence_sha256=digest, original_description=found[0].description))
    assert all(sha(p) == h for p, h in pins.items())
    a.output.mkdir(parents=True, exist_ok=False)
    fasta = a.output / 'proteins.faa'
    fasta.write_text(''.join('>' + r['label'] + '\n' + r['sequence'] + '\n' for r in records))
    for r in records:
        del r['sequence']
    receipt = dict(status='complete_focal_source_sequence_export', source_hashes=pins,
                   proteins=records, artifacts={fasta.name: sha(fasta)},
                   scope='Exact source sequences matching domain annotation hashes. '
                   'Does not validate completeness, alignment, orthology, or domain events.')
    (a.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(records, indent=2))


if __name__ == '__main__':
    main()
