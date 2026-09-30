#!/usr/bin/env python3
"""Plot expanded screening attrition against full event and manifest denominators."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
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
                raise ValueError('Changed attrition figure source: ' + path)

    verify()
    root = Path(plan['source'])
    receipt = json.loads((root / 'receipt.json').read_text())
    audit = json.loads(Path(plan['audit']).read_text())
    assert audit['status'] == 'passed_full_expanded_pair_event_and_taxon_coverage_readback'
    assert audit['producer_receipt_sha256'] == sha(root / 'receipt.json')
    assert audit['pair_mask_rows'] == receipt['pair_mask_rows'] == 269624
    assert audit['event_rows'] == receipt['events'] == 935353
    assert audit['event_mask_rows'] == receipt['event_mask_rows'] == 1870706
    assert audit['taxon_screen_rows'] == receipt['taxon_screen_rows'] == 12624
    assert audit['event_pass_counts'] == receipt['event_pass_counts']
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest
    data = pd.read_csv(root / 'taxon_screen_coverage.tsv', sep='\t', keep_default_na=False)
    keys = ['guide', 'taxon_id', 'mask', 'screen']
    assert len(data) == 12624 and not data.duplicated(keys).any()
    assert data.taxon_id.nunique() == 526
    assert (data.passed_events <= data.distinct_model_pair).all()
    assert (data.distinct_model_pair <= data.both_models).all() and (data.both_models <= data.events).all()
    summary = []
    for (guide, mask, screen), sub in data.groupby(['guide', 'mask', 'screen'], sort=True):
        assert len(sub) == 526 and sub.taxon_id.nunique() == 526
        passed, denominator = int(sub.passed_events.sum()), int(sub.events.sum())
        assert denominator == {'profile': 467663, 'mafft': 467690}[guide]
        assert passed == receipt['event_pass_counts'][guide + ':' + mask + ':' + screen]
        summary.append(dict(guide=guide, mask=mask, screen=screen, candidate_events=denominator,
                            both_models=int(sub.both_models.sum()), distinct_models=int(sub.distinct_model_pair.sum()),
                            passed_events=passed, pass_fraction_of_all_events=passed / denominator,
                            manifest_taxa=len(sub), taxa_with_candidate_events=int(sub.events.gt(0).sum()),
                            taxa_with_both_models=int(sub.both_models.gt(0).sum()),
                            taxa_with_passed_events=int(sub.passed_events.gt(0).sum())))
    frame = pd.DataFrame(summary)
    assert len(frame) == 24
    data['broad_lineage'] = data.lineage.str.split(';').str[0]
    lineage = data.groupby(['guide', 'mask', 'screen', 'broad_lineage'])[
        ['events', 'both_models', 'distinct_model_pair', 'passed_events']].sum()
    lineage['manifest_taxa'] = data.groupby(['guide', 'mask', 'screen', 'broad_lineage']).size()
    lineage['taxa_with_passed_events'] = data[data.passed_events.gt(0)].groupby(
        ['guide', 'mask', 'screen', 'broad_lineage']).size().reindex(lineage.index, fill_value=0)
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    frame.to_csv(out / 'screen_attrition_summary.tsv', sep='\t', index=False)
    lineage.to_csv(out / 'lineage_screen_attrition.tsv', sep='\t')
    pd.testing.assert_frame_equal(frame, pd.read_csv(out / 'screen_attrition_summary.tsv', sep='\t', float_precision='round_trip'))
    pd.testing.assert_frame_equal(lineage, pd.read_csv(out / 'lineage_screen_attrition.tsv', sep='\t').set_index(lineage.index.names))
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9))
    screen_names = [s['id'] for s in receipt['screens']]
    combinations = [('profile', 'full'), ('profile', 'plddt70'), ('mafft', 'full'), ('mafft', 'plddt70')]
    for ix, (guide, mask) in enumerate(combinations):
        sub = frame[frame.guide.eq(guide) & frame['mask'].eq(mask)].set_index('screen').loc[screen_names]
        x = np.arange(len(screen_names)) + (ix - 1.5) * .2
        label = guide.title() + ', ' + ('full protein' if mask == 'full' else 'pLDDT ≥70')
        color = '#0072B2' if guide == 'profile' else '#D55E00'
        hatch = '///' if mask == 'plddt70' else None
        axes[0].bar(x, 100 * sub.pass_fraction_of_all_events, width=.19, color=color,
                    alpha=.75, hatch=hatch, label=label)
        axes[1].bar(x, sub.taxa_with_passed_events, width=.19, color=color, alpha=.75, hatch=hatch)
    labels = [f"{s['minimum_aligned_residues']} aa\n{int(100*s['minimum_original_coverage'])}%" for s in receipt['screens']]
    for ax in axes:
        ax.set_xticks(np.arange(len(labels)), labels)
        ax.set_xlabel('Minimum aligned residues / coverage of each original protein')
    axes[0].set_ylabel('Passing events (% of full terminal candidate ledger)')
    axes[0].set_title('A  Event attrition', loc='left')
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].set_ylabel('Taxa with ≥1 passing event (of 526 entries)')
    axes[1].set_title('B  Taxon representation', loc='left')
    axes[1].set_ylim(0, 526)
    fig.suptitle('Expanded duplication coverage and confidence-mask sensitivity', fontsize=13)
    fig.text(.5, .015, 'Both input orders must pass. Missing/same-model events remain in denominators. Guides and events overlap; descriptive screening, not biological inference.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .065, 1, .95))
    for suffix in ['png', 'pdf']:
        fig.savefig(out / ('expanded_duplication_screen_attrition.' + suffix), dpi=180)
    plt.close(fig)
    verify()
    result = dict(status='complete_source_checked_expanded_coverage_figure_pending_visual_review',
                  plan_sha256=sha(args.plan), source_hashes=pins, summary_rows=24, lineage_rows=len(lineage),
                  event_rows=935353, manifest_taxa=526, artifacts={p.name: sha(p) for p in out.iterdir()},
                  scope='Every screening aggregate checked against complete independently audited taxon rows and full event counts; all six screens, both masks and guides retained. '
                        'Figures use the entire terminal singleton-side candidate ledger and all manifest entries, including missing/same-model and zero-event taxa. '
                        'Broad-lineage labels are descriptive manifest annotations; no confidence calibration, missingness correction, independent-event test or biological effect claim.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
