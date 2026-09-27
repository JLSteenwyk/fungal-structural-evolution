#!/usr/bin/env python3
"""Enumerate incompatible old/new internal splits without support cutoffs."""
import argparse
import csv
import hashlib
import itertools
import json
from collections import defaultdict
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    root = Path('results/cds/codon-alignment-tree-comparison-20260927-v1')
    proof_path = Path('metadata/codon_tree_comparison_completed_readback_20260927.json')
    proof = json.loads(proof_path.read_text())
    receipt = json.loads((root / 'receipt.json').read_text())
    assert proof['status'] == 'passed_full_codon_tree_comparison_independent_readback'
    assert proof['source_receipt_sha256'] == sha(root / 'receipt.json')
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest
    cases = {r['case_id']: r for r in rows(root / 'cases.tsv')}
    groups = defaultdict(list)
    for r in rows(root / 'splits.tsv'):
        groups[r['case_id']].append(r)
    output = []
    candidate_count = 0
    for case, splits in sorted(groups.items()):
        taxa = {r['smaller_side_taxa'] for r in splits if r['internal'] == 'False'}
        assert len(taxa) == int(cases[case]['taxa'])
        old = [r for r in splits if r['internal'] == 'True' and r['status'] == 'original_only']
        new = [r for r in splits if r['internal'] == 'True' and r['status'] == 'local_only']
        for a, b in itertools.product(old, new):
            candidate_count += 1
            x, y = set(a['smaller_side_taxa'].split(';')), set(b['smaller_side_taxa'].split(';'))
            cells = [x & y, x - y, y - x, taxa - (x | y)]
            if not all(cells):
                continue
            output.append(dict(case_id=case, original_split=a['smaller_side_taxa'], local_split=b['smaller_side_taxa'],
                               original_ufb=float(a['original_ufb']), local_ufb=float(b['local_ufb']),
                               minimum_ufb=min(float(a['original_ufb']), float(b['local_ufb'])),
                               original_sh_alrt=float(a['original_sh_alrt']), local_sh_alrt=float(b['local_sh_alrt']),
                               original_length=float(a['original_length']), local_length=float(b['local_length']),
                               witness_quartet=';'.join(min(c) for c in cells),
                               historical_review_flags=cases[case]['historical_review_flags'],
                               alignment_median_pair_retention=cases[case]['alignment_median_pair_retention']))
    assert len({r['case_id'] for r in output}) == proof['changed_topology_cases']
    args.output.mkdir(parents=True, exist_ok=False)
    path = args.output / 'incompatible_splits.tsv'
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(output)
    # An entire empirical support curve avoids a selected significance threshold.
    curve = []
    for support in sorted({r['minimum_ufb'] for r in output}):
        selected = [r for r in output if r['minimum_ufb'] >= support]
        curve.append(dict(minimum_both_ufb=support, incompatible_pairs=len(selected), cases=len({r['case_id'] for r in selected})))
    with (args.output / 'support_curve.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(curve[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(curve)
    report = dict(status='complete_codon_split_conflict_enumeration_pending_readback',
                  source_receipt_sha256=sha(root / 'receipt.json'), source_readback_sha256=sha(proof_path),
                  script_sha256=sha(__file__), candidate_pairs=candidate_count,
                  incompatible_pairs=len(output), cases=len({r['case_id'] for r in output}),
                  zero_length_pairs=sum(r['original_length'] == 0 or r['local_length'] == 0 for r in output),
                  maximum_joint_ufb=max(r['minimum_ufb'] for r in output),
                  artifacts={p.name: sha(p) for p in args.output.iterdir()},
                  scope='Four nonempty split intersections prove unrooted incompatibility. Every old-only/new-only pair is considered. Support curves are descriptive, not significance, independence, biological validity, or selection eligibility. Zero-length splits and historical flags retained.')
    (args.output / 'receipt.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
