#!/usr/bin/env python3
"""Prepare the next length tier of experimental controls without truncation."""
import argparse
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from audit_busco_gene_copies import read_table, sha
from assess_pae_sensitivity import checked_receipt
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--max-length', type=int, required=True)
    a = p.parse_args()
    if a.output.exists() or a.max_length <= 512:
        raise ValueError('Require a fresh output and a longer length tier')
    source = Path('data/prediction_inputs/experimental-controls-v1')
    checked_receipt(source)
    fasta = Path('data/domains/marker-inputs-v1/sequences.faa')
    disposition = read_table(source / 'sequence_disposition.tsv')
    targets = {r['sequence_id']: r for r in disposition if r['status'] == 'deferred_length'}
    sequences = {r.id: str(r.seq) for r in SeqIO.parse(fasta, 'fasta') if r.id in targets}
    if set(sequences) != set(targets):
        raise ValueError('Missing reference sequences')
    reserved, pins = {}, {}
    for f in sorted(Path('data/prediction_inputs').glob('*/candidates.faa')):
        receipt = json.loads((f.parent / 'receipt.json').read_text())
        if sha(f) != receipt['artifacts']['candidates.faa']:
            raise ValueError('Changed reservation FASTA: ' + str(f))
        pins[str(f)] = sha(f)
        for r in SeqIO.parse(f, 'fasta'):
            reserved.setdefault(r.id, []).append(f.parent.name)
    chosen, rows = {}, []
    for sid, sequence in sorted(sequences.items()):
        if sid != 'S' + hashlib.sha256(sequence.encode()).hexdigest() or len(sequence) != int(targets[sid]['length']):
            raise ValueError('Reference identity/length mismatch')
        status = ('reserved_existing_queue' if sid in reserved else
                  'deferred_noncanonical' if not set(sequence) <= set('ACDEFGHIKLMNPQRSTVWY') else
                  'deferred_length' if len(sequence) > a.max_length else 'new_control')
        if status == 'new_control':
            chosen[sid] = sequence
        rows.append(dict(sequence_id=sid, length=len(sequence), status=status,
                         existing_queues=';'.join(reserved.get(sid, []))))
    if not chosen:
        raise ValueError('No unreserved eligible controls')
    a.output.mkdir(parents=True)
    (a.output / 'candidates.faa').write_text(''.join('>' + sid + '\n' + chosen[sid] + '\n'
        for sid in sorted(chosen, key=lambda s: (len(chosen[s]), s))))
    if {r.id: str(r.seq) for r in SeqIO.parse(a.output / 'candidates.faa', 'fasta')} != chosen:
        raise ValueError('Written queue differs')
    links = [r for r in read_table(source / 'experimental_reference_links.tsv') if r['sequence_id'] in chosen]
    if {r['sequence_id'] for r in links} != set(chosen):
        raise ValueError('Missing experimental links')
    write_table(a.output / 'experimental_reference_links.tsv', links)
    write_table(a.output / 'sequence_disposition.tsv', rows)
    receipt = dict(status='complete_prediction_queue_preparation', prediction_candidates=len(chosen),
        candidate_residues=sum(map(len, chosen.values())), max_length=a.max_length,
        original_deferred_sequences=len(targets), reference_links=len(links),
        source_receipt_sha256=sha(source / 'receipt.json'), source_fasta_sha256=sha(fasta),
        reservation_fasta_sha256=pins, script_sha256=sha(Path(__file__)),
        selection='All previously length-deferred experimental controls within this tier, excluding every existing candidate queue; full canonical sequences, no truncation or agreement-based selection.',
        artifacts={f.name: sha(f) for f in a.output.iterdir()})
    (a.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if not isinstance(v, dict)}, indent=2))


if __name__ == '__main__':
    main()
