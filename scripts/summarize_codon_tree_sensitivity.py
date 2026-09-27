#!/usr/bin/env python3
"""Summarize the complete independently verified local/original tree comparison."""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path('results/cds/codon-alignment-tree-comparison-20260927-v1')
    proof_path = Path('metadata/codon_tree_comparison_completed_readback_20260927.json')
    proof = json.loads(proof_path.read_text())
    source = json.loads((root / 'receipt.json').read_text())
    assert proof['status'] == 'passed_full_codon_tree_comparison_independent_readback'
    assert proof['source_receipt_sha256'] == sha(root / 'receipt.json')
    for name, expected in source['artifacts'].items():
        assert sha(root / name) == expected
    cases = pd.read_csv(root / 'cases.tsv', sep='\t', keep_default_na=False)
    splits = pd.read_csv(root / 'splits.tsv', sep='\t')
    matched = cases[cases.comparison_status == 'matched'].copy()
    matched['topology_changed'] = matched.rf_distance.astype(int) > 0
    matched['has_historical_flags'] = matched.historical_review_flags != ''
    changed = matched[matched.topology_changed]
    assert len(changed) == proof['changed_topology_cases']
    assert len(matched) == proof['dispositions']['matched']
    flags = matched.groupby(['topology_changed', 'has_historical_flags']).size().rename('cases').reset_index()
    summaries = []
    support_rows = []
    for status, column in [('original_only', 'original_ufb'), ('local_only', 'local_ufb')]:
        values = splits.loc[splits.internal & (splits.status == status), column]
        for support, count in values.value_counts().sort_index().items():
            support_rows.append(dict(split_status=status, ufb_percent=support, branches=count))
        summaries.append(dict(split_status=status, branches=len(values),
                              ufb_quantiles={str(q): float(v) for q, v in values.quantile([0, .25, .5, .75, 1]).items()}))
    args.output.mkdir(parents=True, exist_ok=False)
    flags.to_csv(args.output / 'topology_by_historical_flags.tsv', sep='\t', index=False)
    pd.DataFrame(support_rows).to_csv(args.output / 'changed_split_support_distribution.tsv', sep='\t', index=False)
    changed.to_csv(args.output / 'changed_cases.tsv', sep='\t', index=False)
    receipt = dict(status='complete_verified_codon_tree_sensitivity_summary',
                   comparison_receipt_sha256=sha(root / 'receipt.json'), comparison_readback_sha256=sha(proof_path),
                   script_sha256=sha(__file__), ledger_cases=len(cases), matched_cases=len(matched),
                   changed_cases=len(changed), changed_fraction=len(changed)/len(matched),
                   changed_with_historical_flags=int(changed.has_historical_flags.sum()),
                   changed_normalized_rf_quantiles={str(q): float(v) for q, v in changed.normalized_rf.astype(float).quantile([0, .25, .5, .75, 1]).items()},
                   changed_split_support=summaries,
                   case_median_retention_by_topology_change={str(k): float(v) for k, v in matched.assign(retention=matched.alignment_median_pair_retention.astype(float)).groupby('topology_changed').retention.median().items()},
                   artifacts={p.name: sha(p) for p in args.output.iterdir()},
                   scope='Descriptive full matched-case sensitivity; all original 1712 ledger dispositions retained in source. Zero-length resolutions retained. Support quantiles are not significance tests; changed cases and branches are correlated. High support on some lost/new splits does not alone establish strongly supported incompatible splits. Historical flags are review triggers, not eligibility decisions.')
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
