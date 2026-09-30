#!/usr/bin/env python3
"""Independently merge full conflict tables and record branch-level method sensitivity."""
import argparse
import json
from pathlib import Path

import pandas as pd
from screen_duplication_alignment_reuse import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    pins = {str(args.plan): sha(args.plan), **plan['pins']}

    def verify():
        for path, digest in pins.items():
            if sha(path) != digest:
                raise ValueError('Changed support-sensitivity source: ' + path)

    verify()
    root = Path(plan['source'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['status'] == 'complete_full_alignment_supported_conflict_sensitivity_with_exhaustive_readbacks'
    assert receipt['paired_rows_checked'] == 261500
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest
    execution_plan = json.loads(Path(plan['source_plan']).read_text())
    assert receipt['plan_sha256'] == sha(plan['source_plan'])
    keys = ['marker', 'guide', 'full_split_taxa', 'sh_alrt_cutoff']
    frames = []
    for label in ['profile', 'mafft']:
        spec = execution_plan['sources'][label]
        proof = json.loads(Path(spec['proof']).read_text())
        diagnostic = Path(spec['diagnostic'])
        dr = json.loads((diagnostic / 'receipt.json').read_text())
        assert proof['status'] == 'passed_full_exhaustive_supported_marker_conflict_readback'
        assert proof['diagnostic_receipt_sha256'] == sha(diagnostic / 'receipt.json')
        assert proof['grid_status_and_witness_rows_checked'] == 261500
        assert proof['independent_full_set_conflict_maxima_checked'] == 130750
        for name, digest in dr['artifacts'].items():
            assert sha(diagnostic / name) == digest
        frame = pd.read_csv(diagnostic / 'marker_guide_conflict.tsv', sep='\t', keep_default_na=False)
        assert len(frame) == 261500 and not frame.duplicated(keys).any()
        frame = frame[keys + ['marker_taxa', 'status']].rename(columns={
            'marker_taxa': label + '_marker_taxa', 'status': label + '_status'})
        frames.append(frame)
    expected = frames[0].merge(frames[1], on=keys, how='outer', validate='one_to_one', indicator=True)
    assert len(expected) == 261500 and expected['_merge'].eq('both').all()
    expected = expected.drop(columns='_merge')
    expected['classification_changed'] = expected.profile_status.ne(expected.mafft_status).astype(int)
    actual = pd.read_csv(root / 'paired_guide_conflict.tsv', sep='\t', keep_default_na=False)
    assert not actual.duplicated(keys).any()
    pd.testing.assert_frame_equal(actual[expected.columns].set_index(keys).sort_index(),
                                  expected.set_index(keys).sort_index())
    assert int(expected.classification_changed.sum()) == receipt['classification_changed_rows']
    transition_keys = ['guide', 'sh_alrt_cutoff', 'profile_status', 'mafft_status']
    counts = expected.groupby(transition_keys).size().rename('marker_edge_cells').reset_index()
    stored = pd.read_csv(root / 'classification_transition_counts.tsv', sep='\t')
    pd.testing.assert_frame_equal(counts.set_index(transition_keys).sort_index(),
                                  stored.set_index(transition_keys).sort_index())
    assert len(counts) == receipt['transition_rows'] == 52
    membership = pd.read_csv(plan['membership'], sep='\t', keep_default_na=False)
    assert len(membership) == 125 and membership.marker.is_unique
    membership['same_original_tips'] = membership.first_only_taxa.eq('') & membership.second_only_taxa.eq('')
    assert int(membership.same_original_tips.sum()) == 95
    expected = expected.merge(membership[['marker', 'same_original_tips']], on='marker',
                              how='left', validate='many_to_one')
    assert expected.same_original_tips.notna().all()
    expected['stable_concordance'] = expected.profile_status.eq('supported_concordance') & expected.mafft_status.eq('supported_concordance')
    expected['stable_conflict'] = expected.profile_status.eq('supported_conflict') & expected.mafft_status.eq('supported_conflict')
    expected['concordance_conflict_reversal'] = (
        expected.profile_status.eq('supported_concordance') & expected.mafft_status.eq('supported_conflict')) | (
        expected.profile_status.eq('supported_conflict') & expected.mafft_status.eq('supported_concordance'))
    metrics = ['classification_changed', 'stable_concordance', 'stable_conflict', 'concordance_conflict_reversal']
    edge_keys = ['guide', 'full_split_taxa', 'sh_alrt_cutoff']
    edge = expected.groupby(edge_keys)[metrics].sum()
    edge['marker_edge_cells'] = expected.groupby(edge_keys).size()
    assert len(edge) == 2092 and edge.marker_edge_cells.eq(125).all()
    strata_keys = ['guide', 'sh_alrt_cutoff', 'same_original_tips']
    strata = expected.groupby(strata_keys)[metrics].sum()
    strata['marker_edge_cells'] = expected.groupby(strata_keys).size()
    assert len(strata) == 8
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    for name, frame in [('branch_method_sensitivity.tsv', edge),
                        ('membership_stratified_sensitivity.tsv', strata)]:
        frame.to_csv(out / name, sep='\t')
        pd.testing.assert_frame_equal(pd.read_csv(out / name, sep='\t').set_index(frame.index.names), frame)
    verify()
    result = dict(status='passed_full_supported_conflict_serialized_join_and_branch_summary_readback',
                  plan_sha256=sha(args.plan), source_receipt_sha256=sha(root / 'receipt.json'),
                  source_hashes=pins, paired_rows_checked=261500, transition_rows_checked=52,
                  branch_summary_rows=2092, membership_strata=8, same_tip_markers=95, different_tip_markers=30,
                  classification_changed_rows=int(expected.classification_changed.sum()),
                  artifacts={p.name: sha(p) for p in out.iterdir()}, scientific_eligibility=False,
                  scope='Independent pandas outer one-to-one merge reconstructs every serialized paired classification and all transition counts from two exhaustively audited source tables. '
                        'Branch and original-tip-membership summaries retain all cells and cutoffs; summary serialization checked. '
                        'No native inference rerun, sampling independence, alignment-only causal interpretation, species-tree acceptance or significance claim.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
