#!/usr/bin/env python3
"""Screen within-cluster Pfam pairs for shared-model and annotation support."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
from itertools import combinations, groupby
import json
from pathlib import Path
import re
import sqlite3


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8388608), b''):
            h.update(block)
    return h.hexdigest()


def pair_counts(left, right):
    shared = len(left & right)
    return (len(left), len(right), shared, len(left - right), len(right - left),
            len(left) * len(right) - shared,
            len(left - right) * len(right - left))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    plan_hash = sha(args.plan)

    def verify():
        if sha(args.plan) != plan_hash:
            raise ValueError('Plan changed')
        for path, digest in plan['pins'].items():
            if sha(path) != digest:
                raise ValueError('Input changed: ' + path)

    verify()
    receipt = json.loads(Path(plan['source_receipt']).read_text())
    if receipt['status'] != 'complete_domain_cluster_pfam_source_set_sql_agreement':
        raise ValueError('Completed annotation join required')
    if receipt['artifacts']['cluster_pfam.sqlite'] != plan['pins'][plan['database']]:
        raise ValueError('Database/receipt binding differs')
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    db = sqlite3.connect('file:' + plan['database'] + '?mode=ro', uri=True)
    counts = Counter(); represented = Counter(); multi = Counter(); rows_seen = 0
    table = out / 'pfam_pair_screen.tsv'
    header = ['boundary', 'representative', 'pfam_left', 'pfam_right',
              'clan_left', 'clan_right', 'clan_relation', 'models_left',
              'models_right', 'models_shared', 'models_left_only',
              'models_right_only', 'distinct_model_pairs',
              'exclusive_model_pairs', 'shared_intervals']
    with table.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        query = ('SELECT boundary,representative,pfam_accession,pfam_clan,model,interval_id '
                 'FROM links ORDER BY boundary,representative,pfam_accession,pfam_clan')
        for (boundary, rep), group in groupby(db.execute(query), lambda r: r[:2]):
            represented[boundary] += 1
            models = defaultdict(set); intervals = defaultdict(set); clans = {}
            for _, _, accession, clan, model, interval in group:
                if clan and not re.fullmatch(r'CL\d{4}', clan):
                    raise ValueError('Unexpected clan identity')
                if accession in clans and clans[accession] != clan:
                    raise ValueError('Inconsistent accession/clan mapping')
                clans[accession] = clan
                models[accession].add(model); intervals[accession].add(interval)
                rows_seen += 1
            if len(models) < 2:
                continue
            multi[boundary] += 1
            for left, right in combinations(sorted(models), 2):
                cl, cr = clans[left], clans[right]
                relation = ('unknown_clan' if not cl or not cr else
                            'same_known_clan' if cl == cr else 'different_known_clans')
                values = pair_counts(models[left], models[right])
                overlap = len(intervals[left] & intervals[right])
                writer.writerow([boundary, rep, left, right, cl, cr, relation,
                                 *values, overlap])
                counts[boundary + ':pairs'] += 1
                counts[boundary + ':' + relation] += 1
                if values[5]: counts[boundary + ':distinct_model_support'] += 1
                if values[6]: counts[boundary + ':exclusive_model_support'] += 1
                if overlap: counts[boundary + ':shared_interval'] += 1
    if rows_seen != receipt['boundary_links']:
        raise ValueError('Incomplete source links')
    if dict(represented) != receipt['represented_clusters_by_boundary']:
        raise ValueError('Incomplete represented clusters')
    if dict(multi) != receipt['multiple_pfam_clusters_by_boundary']:
        raise ValueError('Incomplete multi-Pfam clusters')
    # Read every output row and independently enumerate the relevant model pairs.
    read = 0
    with table.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            sets = []
            for accession in [row['pfam_left'], row['pfam_right']]:
                sets.append({r[0] for r in db.execute(
                    'SELECT DISTINCT model FROM links INDEXED BY group_idx WHERE boundary=? AND representative=? AND pfam_accession=?',
                    (row['boundary'], row['representative'], accession))})
            left, right = sets
            # Enumeration avoids reusing the production count formula.
            distinct = sum(a != b for a in left for b in right)
            exclusive = sum(a not in right and b not in left for a in left for b in right)
            if distinct != int(row['distinct_model_pairs']) or exclusive != int(row['exclusive_model_pairs']):
                raise ValueError('Model-pair enumeration differs')
            read += 1
    if read != sum(v for k, v in counts.items() if k.endswith(':pairs')):
        raise ValueError('Serialized row count differs')
    db.close(); verify()
    result = dict(status='complete_domain_pfam_pair_screen', plan_sha256=plan_hash,
                  source_links=rows_seen, multi_pfam_clusters=dict(multi),
                  pair_rows=read, counts=dict(sorted(counts.items())),
                  artifacts={table.name: sha(table)},
                  scope='All within-cluster accession pairs, separated by boundary view. '
                  'Distinct model pairs exclude identical source models; exclusive model pairs '
                  'require each model to carry only its own accession among this pair within this cluster. '
                  'Shared intervals detect identical intervals, not all overlaps. '
                  'No confidence qualification, independent clustering, direct structural alignment, '
                  'homology, functional equivalence, or evolutionary inference is established.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
