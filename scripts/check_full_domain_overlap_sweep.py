#!/usr/bin/env python3
"""Compare the overlap sweep to exhaustive integer-residue intersections."""
import argparse
import itertools
import json
from pathlib import Path
import random
from map_full_domain_overlaps import pairs, sha


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    rng = random.Random(20260917)
    nested = {('PF_A', 'PF_B'), ('PF_C', 'PF_A')}
    checked = 0
    for case in range(400):
        hits = []
        for i in range(case % 36):
            start = rng.randint(1, 100)
            end = start + rng.randint(0, 50)
            left = rng.randint(start, end)
            right = rng.randint(left, end)
            hits.append(dict(hit_id=f'H{i}', pfam_accession=rng.choice(['PF_A','PF_B','PF_C']),
                             pfam_clan=rng.choice(['','CL_A','CL_B']), envelope_start=start,
                             envelope_end=end, alignment_start=left, alignment_end=right))
        expected = {}
        for x, y in itertools.combinations(hits, 2):
            xe = set(range(x['envelope_start'], x['envelope_end'] + 1))
            ye = set(range(y['envelope_start'], y['envelope_end'] + 1))
            xa = set(range(x['alignment_start'], x['alignment_end'] + 1))
            ya = set(range(y['alignment_start'], y['alignment_end'] + 1))
            if xe & ye:
                key = frozenset([x['hit_id'],y['hit_id']])
                expected[key] = (len(xe & ye), len(xa & ya),
                                 int(x['pfam_accession'] == y['pfam_accession']),
                                 int(bool(x['pfam_clan']) and x['pfam_clan'] == y['pfam_clan']),
                                 int(((x['pfam_accession'], y['pfam_accession']) in nested and ya <= xa)
                                     or ((y['pfam_accession'], x['pfam_accession']) in nested and xa <= ya)))
        rng.shuffle(hits)
        observed = {}
        for row in pairs(hits, nested):
            key = frozenset([row['hit_a'],row['hit_b']])
            assert key not in observed
            observed[key] = tuple(row[k] for k in ['envelope_overlap','alignment_overlap','same_family','same_clan','curated_and_coordinate_nesting'])
        assert observed == expected, case
        checked += len(expected)
    receipt = dict(status='passed_exhaustive_residue_set_overlap_fixtures', cases=400,
                   overlapping_pairs_checked=checked, script_sha256=sha(Path(__file__)),
                   sweep_sha256=sha(Path('scripts/map_full_domain_overlaps.py')),
                   scope='Exact envelope/alignment overlap sizes, unique pair universe, clan/family flags and directed model-plus-coordinate nesting checked against residue sets under shuffled inputs. Software fixtures, not biological validation.')
    a.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
