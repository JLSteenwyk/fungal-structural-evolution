#!/usr/bin/env python3
"""Freeze full-marker prediction candidates while distinguishing unfinished reuse searches."""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict, deque
from pathlib import Path
from Bio import SeqIO
from prepare_pfam import ROOT, digest


def complete_json_lines(path):
    data = path.read_bytes()
    # Concurrent appenders may have one unfinished last line. Freeze complete records only.
    return data[:data.rfind(b'\n') + 1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use an immutable new queue snapshot')
    inputs = ROOT / 'data/domains/marker-inputs-v1'
    ir = json.loads((inputs / 'receipt.json').read_text())
    for name, expected in ir['artifacts'].items():
        if digest(inputs / name) != expected:
            raise ValueError('Changed marker input')
    sequences = {r.id: str(r.seq) for r in SeqIO.parse(inputs / 'sequences.faa', 'fasta')}
    links = list(csv.DictReader((inputs / 'protein_links.tsv').open(), delimiter='\t'))
    taxa = defaultdict(set)
    for row in links:
        sid = row['sequence_id']
        if hashlib.sha256(sequences[sid].encode()).hexdigest() != row['sequence_sha256']:
            raise ValueError('Marker sequence identity mismatch')
        taxa[sid].add(row['taxon_id'])
    match_bytes = complete_json_lines(ROOT / 'data/raw/uniprot_matches.jsonl')
    matches = {}
    for line in match_bytes.splitlines():
        row = json.loads(line)
        if row['status'] == 'matched':
            matches[row['taxon_id']] = row
    nominated, match_artifacts = set(), {}
    for row in matches.values():
        path = ROOT / row['match_path']
        before = digest(path)
        with path.open() as handle:
            nominated.update('S' + r['sequence_sha256'] for r in csv.DictReader(handle, delimiter='\t'))
        if digest(path) != before:
            raise ValueError('UniProt match file changed during snapshot')
        match_artifacts[row['match_path']] = before
    model_bytes = complete_json_lines(ROOT / 'data/raw/afdb_models.jsonl')
    verified = set()
    for line in model_bytes.splitlines():
        row = json.loads(line)
        if row['status'] == 'verified':
            verified.update('S' + m['sequence_sha256'] for m in row['models'])
    statuses, buckets = {}, defaultdict(list)
    for sid, sequence in sequences.items():
        completed_taxa = sorted(taxa[sid] & matches.keys())
        status = ('verified_reuse_receipt' if sid in verified else 'nominated_reuse_pending'
                  if sid in nominated else 'no_candidate_in_completed_inventory'
                  if completed_taxa else 'inventory_pending')
        statuses[sid] = status
        if status == 'no_candidate_in_completed_inventory':
            buckets[completed_taxa[0]].append(sid)
    queues = {taxon: deque(sorted(ids, key=lambda s: (len(sequences[s]), s)))
              for taxon, ids in buckets.items()}
    ordered = []
    while any(queues.values()):
        for taxon in sorted(queues):
            if queues[taxon]:
                ordered.append(queues[taxon].popleft())
    args.output.mkdir(parents=True)
    with (args.output / 'candidates.faa').open('w') as handle:
        for sid in ordered:
            handle.write('>' + sid + '\n' + sequences[sid] + '\n')
    with (args.output / 'all_marker_links.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, list(links[0]) + ['reuse_state'], delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(row | {'reuse_state': statuses[row['sequence_id']]} for row in links)
    (args.output / 'uniprot_snapshot.jsonl').write_bytes(match_bytes)
    (args.output / 'afdb_snapshot.jsonl').write_bytes(model_bytes)
    receipt = {'status': 'complete_prediction_queue_preparation',
        'marker_input_receipt_sha256': digest(inputs / 'receipt.json'),
        'script_sha256': digest(Path(__file__)), 'completed_inventory_taxa': len(matches),
        'all_marker_sequences': len(sequences), 'all_marker_links': len(links),
        'unique_sequence_states': dict(Counter(statuses.values())),
        'prediction_candidates': len(ordered),
        'candidate_residues': sum(len(sequences[s]) for s in ordered),
        'candidates_at_most_512_residues': sum(len(sequences[s]) <= 512 for s in ordered),
        'candidate_taxa': len(set().union(*(taxa[s] for s in ordered))) if ordered else 0,
        'match_artifacts': match_artifacts,
        'ordering': 'Round-robin over lexically first completed taxon for each exact sequence, shortest first within each taxon.',
        'interpretation': 'No candidate in the completed taxon-specific UniProt cross-reference inventory is not proof of absence from all structural databases. Unfinished nominated retrievals and unfinished inventories are distinct. All taxon links retained; later inventories require a new snapshot.',
        'artifacts': {p.name: digest(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'match_artifacts'}, indent=2))


if __name__ == '__main__':
    main()
