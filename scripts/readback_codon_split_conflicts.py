#!/usr/bin/env python3
"""Check every candidate conflict by enumerating discordantly resolved quartets."""
import csv
import hashlib
import itertools
import json
from collections import defaultdict
from pathlib import Path


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    with Path(p).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def resolution(q, side):
    a = tuple(t for t in q if t in side)
    b = tuple(t for t in q if t not in side)
    return tuple(sorted((a, b))) if len(a) == len(b) == 2 else None


def main():
    source = Path('results/cds/codon-alignment-tree-comparison-20260927-v1')
    target = Path('results/cds/codon-split-conflicts-20260927-v1')
    receipt = json.loads((target / 'receipt.json').read_text())
    assert receipt['source_receipt_sha256'] == sha(source / 'receipt.json')
    for root in (source, target):
        for name, digest in json.loads((root / 'receipt.json').read_text())['artifacts'].items():
            assert sha(root / name) == digest
    groups = defaultdict(list)
    for r in read(source / 'splits.tsv'):
        groups[r['case_id']].append(r)
    cases = {r['case_id']: r for r in read(source / 'cases.tsv')}
    actual_rows = read(target / 'incompatible_splits.tsv')
    actual = {(r['case_id'], r['original_split'], r['local_split']): r for r in actual_rows}
    assert len(actual) == len(actual_rows)
    expected = set()
    candidates = 0
    for case, splits in groups.items():
        taxa = sorted(r['smaller_side_taxa'] for r in splits if r['internal'] == 'False')
        for a in (r for r in splits if r['status'] == 'original_only'):
            for b in (r for r in splits if r['status'] == 'local_only'):
                candidates += 1
                x, y = set(a['smaller_side_taxa'].split(';')), set(b['smaller_side_taxa'].split(';'))
                discordant = False
                for quartet in itertools.combinations(taxa, 4):
                    old, new = resolution(quartet, x), resolution(quartet, y)
                    if old is not None and new is not None and old != new:
                        discordant = True
                        break
                key = (case, a['smaller_side_taxa'], b['smaller_side_taxa'])
                if not discordant:
                    assert key not in actual
                    continue
                expected.add(key)
                r = actual[key]
                q = tuple(sorted(r['witness_quartet'].split(';')))
                assert len(set(q)) == 4 and set(q) <= set(taxa)
                assert resolution(q, x) and resolution(q, y) and resolution(q, x) != resolution(q, y)
                for name in ('original_ufb', 'original_sh_alrt', 'original_length'):
                    assert float(r[name]) == float(a[name])
                for name in ('local_ufb', 'local_sh_alrt', 'local_length'):
                    assert float(r[name]) == float(b[name])
                assert float(r['minimum_ufb']) == min(float(a['original_ufb']), float(b['local_ufb']))
                for name in ('historical_review_flags', 'alignment_median_pair_retention'):
                    assert r[name] == cases[case][name]
    assert expected == set(actual)
    curve = read(target / 'support_curve.tsv')
    assert [float(r['minimum_both_ufb']) for r in curve] == sorted({float(r['minimum_ufb']) for r in actual_rows})
    for r in curve:
        selected = [v for v in actual_rows if float(v['minimum_ufb']) >= float(r['minimum_both_ufb'])]
        assert int(r['incompatible_pairs']) == len(selected)
        assert int(r['cases']) == len({v['case_id'] for v in selected})
    assert receipt['candidate_pairs'] == candidates
    assert receipt['incompatible_pairs'] == len(expected)
    assert receipt['cases'] == len({k[0] for k in expected})
    assert receipt['maximum_joint_ufb'] == max(float(v['minimum_ufb']) for v in actual_rows)
    assert receipt['zero_length_pairs'] == sum(float(v['original_length']) == 0 or float(v['local_length']) == 0 for v in actual_rows)
    proof = dict(status='passed_full_quartet_codon_split_conflict_readback',
                 receipt_sha256=sha(target / 'receipt.json'), script_sha256=sha(__file__),
                 candidate_pairs=candidates, incompatible_pairs=len(expected), cases=receipt['cases'],
                 maximum_joint_ufb=receipt['maximum_joint_ufb'], zero_length_pairs=receipt['zero_length_pairs'],
                 scope='Every old-only/new-only pair tested for incompatible resolved quartets, independently of producer intersection algorithm. Complete row grid, witnesses, supports, lengths, flags, retention and empirical curve verified against audited split tables. Does not rerun source tree inference.')
    Path('metadata/codon_split_conflicts_readback_20260927.json').write_text(json.dumps(proof, indent=2)+'\n')
    print(json.dumps(proof, indent=2))


if __name__ == '__main__':
    main()
