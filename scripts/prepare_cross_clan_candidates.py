#!/usr/bin/env python3
"""Inventory every different-known-clan pair from the completed domain screen."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sqlite3


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--screen', type=Path, required=True)
    ap.add_argument('--annotations', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    screen_receipt = args.screen / 'receipt.json'
    annotation_receipt = args.annotations / 'receipt.json'
    table = args.screen / 'pfam_pair_screen.tsv'
    database = args.annotations / 'cluster_pfam.sqlite'
    sources = {str(p): sha(p) for p in [screen_receipt, annotation_receipt, table, database, Path(__file__)]}
    sr = json.loads(screen_receipt.read_text()); ar = json.loads(annotation_receipt.read_text())
    if sr['status'] != 'complete_domain_pfam_pair_screen' or ar['status'] != 'complete_domain_cluster_pfam_source_set_sql_agreement':
        raise ValueError('Completed source analyses required')
    if sr['artifacts'][table.name] != sources[str(table)] or ar['artifacts'][database.name] != sources[str(database)]:
        raise ValueError('Artifact binding differs')
    with table.open() as f:
        reader = csv.DictReader(f, delimiter='\t'); header = reader.fieldnames
        selected = [r for r in reader if r['clan_relation'] == 'different_known_clans']
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / 'candidates.tsv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(selected)
    db = sqlite3.connect('file:' + str(database) + '?mode=ro', uri=True)
    models = set(); intervals = set(); accession_pairs = set(); counts = Counter()
    model_pairs = Counter(); member_count = 0
    with (args.output / 'candidate_members.tsv').open('w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['boundary', 'representative', 'pfam_left', 'pfam_right',
                         'side', 'pfam_accession', 'model', 'interval_id',
                         'model_exclusive_within_pair_cluster'])
        for r in selected:
            key = (r['boundary'], r['representative'], r['pfam_left'], r['pfam_right'])
            members = []
            for accession in key[2:]:
                members.append(set(db.execute(
                    'SELECT model,interval_id FROM links INDEXED BY group_idx '
                    'WHERE boundary=? AND representative=? AND pfam_accession=?',
                    (*key[:2], accession))))
            model_sets = [{m for m, _ in entries} for entries in members]
            exclusive = len(model_sets[0] - model_sets[1]) * len(model_sets[1] - model_sets[0])
            if exclusive != int(r['exclusive_model_pairs']):
                raise ValueError('Candidate identity/source model counts differ')
            for side in range(2):
                for model, interval in sorted(members[side]):
                    writer.writerow([*key, ['left', 'right'][side], key[2 + side],
                                     model, interval, int(model not in model_sets[1 - side])])
                    models.add(model); intervals.add(interval); member_count += 1
            counts[key[0]] += 1; model_pairs[key[0]] += exclusive
            accession_pairs.add(key[2:])
    db.close()
    for boundary in ['alignment', 'envelope']:
        if counts[boundary] != sr['counts'][boundary + ':different_known_clans']:
            raise ValueError('Incomplete cross-clan scope')
    for p, digest in sources.items():
        if sha(p) != digest:
            raise ValueError('Source changed')
    receipt = dict(status='complete_cross_clan_candidate_inventory',
                   source_hashes=sources, candidates_by_boundary=dict(counts),
                   unique_accession_pairs=len(accession_pairs), source_models=len(models),
                   intervals=len(intervals), candidate_member_rows=member_count,
                   exclusive_model_pairs_by_boundary=dict(model_pairs),
                   artifacts={p.name: sha(p) for p in args.output.glob('*.tsv')},
                   scope='Complete different-known-clan candidate inventory from one joint partition. '
                   'Boundary views overlap. Model exclusivity is local to a candidate pair and cluster. '
                   'No confidence qualification, direct alignment, homology, functional or evolutionary claims.')
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
