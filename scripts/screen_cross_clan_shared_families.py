#!/usr/bin/env python3
"""Join exclusive cross-clan pairs within assigned families to direct comparisons."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from map_cross_clan_families import sha


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mapping', type=Path, required=True)
    p.add_argument('--comparisons', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    members = a.mapping / 'candidate_proteins.tsv'
    comparisons = a.comparisons / 'bidirectional_pairs.tsv'
    sources = [members, comparisons, a.mapping / 'receipt.json', a.comparisons / 'receipt.json', Path(__file__), Path(__file__).with_name('map_cross_clan_families.py')]
    pins = {str(x): sha(x) for x in sources}
    for table, receipt in [(members, sources[2]), (comparisons, sources[3])]:
        assert json.loads(receipt.read_text())['artifacts'][table.name] == pins[str(table)]
    metrics = {(r['interval_a'], r['interval_b']): r for r in read(comparisons)}
    groups = defaultdict(lambda: defaultdict(list))
    keys = ['boundary', 'representative', 'pfam_left', 'pfam_right']
    for r in read(members):
        if r['model_exclusive_within_pair_cluster'] != '1':
            continue
        for guide in ['profile', 'mafft']:
            key = tuple(r[k] for k in keys) + (guide, r[guide + '_family'])
            groups[key][r['side']].append(r)
    output = []
    for key, sides in sorted(groups.items()):
        for l in sides['left']:
            for r in sides['right']:
                assert l['model'] != r['model']
                pair = tuple(sorted([l['interval_id'], r['interval_id']]))
                m = metrics[pair]
                row = dict(zip(keys + ['guide', 'family'], key))
                for side, member in [('left', l), ('right', r)]:
                    for field in ['interval_id', 'model', 'native_gene_id', 'taxon_id', 'protein_id']:
                        row[field + '_' + side] = member[field]
                row.update(m)
                for threshold in [0, 0.5, 0.8, 0.9]:
                    row['screen_' + str(threshold)] = int(float(m['min_tm']) >= .5 and float(m['min_coverage']) >= .8 and min(float(m['min_joint_confidence']), float(m['min_domain_confidence'])) >= threshold)
                output.append(row)
    assert output, 'No shared-family pairs; add explicit empty-table handling before reuse on that scope'
    a.output.mkdir(parents=True, exist_ok=False)
    target = a.output / 'shared_family_pairs.tsv'
    with target.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(output[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(output)
    for path, digest in pins.items():
        assert sha(path) == digest
    receipt = dict(status='complete_shared_family_direct_comparison_join', source_hashes=pins,
                   rows=len(output), unique_interval_pairs=len({(r['interval_a'], r['interval_b']) for r in output}),
                   artifacts={target.name: sha(target)},
                   scope='Complete exclusive-model within-assigned-family comparisons under both guides. Repeated proteins, guides and boundary views are not independent observations. Descriptive TM/coverage/confidence screens do not establish homology, orthology, transitions or significance; no PAE qualification.')
    (a.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
