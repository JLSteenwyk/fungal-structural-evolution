#!/usr/bin/env python3
"""Rejoin the fixed terminal-event universe to a pinned expanded model catalog."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    plan_hash = sha(args.plan)

    def verify():
        assert sha(args.plan) == plan_hash, 'Plan changed'
        for path, expected in plan['pins'].items():
            assert sha(path) == expected, path

    verify()
    models = {}
    with open(plan['links']) as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            gene = row['taxon_id'] + '_' + row['protein_id']
            assert gene not in models, gene
            models[gene] = row['model_id']
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    totals = defaultdict(Counter)
    taxa = defaultdict(set)
    families = defaultdict(set)
    nrows = 0
    labels = ['neither_model', 'one_model', 'both_models']
    with gzip.open(plan['events'], 'rt') as source, gzip.open(out / 'event_coverage.tsv.gz', 'wt') as dest:
        reader = csv.DictReader(source, delimiter='\t')
        fields = reader.fieldnames + ['new_model_a', 'new_model_b', 'new_coverage_class', 'new_same_model']
        writer = csv.DictWriter(dest, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in reader:
            guide = row['guide']
            assert all(row[g].startswith(row['taxon_id'] + '_') for g in ['gene_a', 'gene_b'])
            old_n = int(bool(row['model_a'])) + int(bool(row['model_b']))
            assert row['coverage_class'] == labels[old_n]
            assert int(row['same_model']) == int(bool(row['model_a']) and row['model_a'] == row['model_b'])
            a, b = models.get(row['gene_a'], ''), models.get(row['gene_b'], '')
            new_n = int(bool(a)) + int(bool(b))
            same = int(bool(a) and a == b)
            row.update(new_model_a=a, new_model_b=b, new_coverage_class=labels[new_n], new_same_model=same)
            writer.writerow(row)
            counts = totals[guide]
            counts['events'] += 1
            counts['old_' + labels[old_n]] += 1
            counts['new_' + labels[new_n]] += 1
            counts['new_same_model'] += same
            counts['newly_both'] += int(old_n < 2 and new_n == 2)
            counts['lost_both'] += int(old_n == 2 and new_n < 2)
            if new_n == 2:
                taxa[guide].add(row['taxon_id'])
                families[guide].add(row['family'])
            nrows += 1
    old = json.loads(Path(plan['event_receipt']).read_text())
    for guide in old['guides']:
        counts = totals[guide['guide']]
        assert counts['events'] == guide['counts']['terminal_singleton_events']
        for label in labels:
            assert counts['old_' + label] == guide['counts'][label]
    assert set(totals) == {g['guide'] for g in old['guides']}
    verify()
    receipt = dict(status='complete_expanded_catalog_terminal_event_join_pending_independent_readback',
                   plan_sha256=plan_hash, event_rows=nrows, catalog_links=len(models),
                   guides={g: dict(c, taxa_with_both=len(taxa[g]), families_with_both=len(families[g])) for g, c in totals.items()},
                   artifacts={'event_coverage.tsv.gz': sha(out / 'event_coverage.tsv.gz')},
                   scope='Fixed previously audited terminal singleton-side event universe joined to exact-sequence catalog. Model availability only; no new confidence qualification, event validation, structural divergence test or ascertainment correction. Profiles and MAFFT are dependent sensitivity guides. Missing is absent from this frozen catalog, not biological absence.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
