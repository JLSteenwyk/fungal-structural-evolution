#!/usr/bin/env python3
"""Plot descriptive source availability for all selected representatives."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('completion', 'png', 'pdf', 'receipt'):
        p.add_argument('--' + name, type=Path, required=True)
    args = p.parse_args()
    assert not any(q.exists() for q in [args.png, args.pdf, args.receipt])
    completion = json.loads(args.completion.read_text())
    assert completion['status'] == 'complete_verified_full_coding_structure_source_coupling'
    assert all(t['actual_terminal_exit_code'] == 0 for t in completion['original_transports'])
    table = Path(completion['taxon_table'])
    assert sha(table) == completion['source_hashes'][str(table)]
    with table.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert len(rows) == len({r['taxon_id'] for r in rows}) == 526
    assert Counter(r['study_role'] for r in rows) == dict(ingroup=501, outgroup=25)
    statuses = [
        'original_target_genome_agrees_and_inherited_strict_translation_exact',
        'original_target_review_or_no_strict_translation_agreement',
        'no_original_target_with_separate_derived_translation_evidence',
        'no_original_target_or_exact_derived_translation_evidence']
    states = ['afdb_only', 'esmfold_only', 'both', 'neither', 'alternative_not_assigned']
    grouped = {}
    recounted = Counter()
    for role in ['ingroup', 'outgroup']:
        counts = Counter()
        for row in rows:
            if row['study_role'] != role:
                continue
            cross = {s + ':' + a: int(row[s + ':' + a]) for s in statuses for a in states}
            assert all(v >= 0 for v in cross.values())
            assert sum(cross.values()) == int(row['source_products'])
            assert sum(v for k, v in cross.items() if not k.endswith(':alternative_not_assigned')) == int(row['selected_representatives'])
            recounted.update({k: v for k, v in cross.items() if v})
            for index, status in enumerate(statuses):
                modeled = sum(cross[status + ':' + a] for a in states[:3])
                assert index != 3 or modeled == 0
                if index < 3:
                    counts[status] += modeled
            counts['no_available_model'] += sum(cross[s + ':neither'] for s in statuses)
            counts['representatives'] += int(row['selected_representatives'])
            counts['entries'] += 1
        assert sum(counts[s] for s in statuses[:3]) + counts['no_available_model'] == counts['representatives']
        grouped[role] = dict(counts)
    assert dict(recounted) == completion['coding_structure_cross_counts']
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    categories = [*statuses[:3], 'no_available_model']
    labels = ['Model + original CDS agreement', 'Model + original CDS review',
              'Model + separate derived evidence', 'No available representative model']
    colors = ['#247773', '#d99230', '#8666a9', '#c7ccd2']
    left = [0.0, 0.0]
    roles = ['ingroup', 'outgroup']
    for category, label, color in zip(categories, labels, colors):
        widths = [100 * grouped[r].get(category, 0) / grouped[r]['representatives'] for r in roles]
        ax.barh([0, 1], widths, left=left, label=label, color=color, height=0.48)
        for i, width in enumerate(widths):
            if width >= 8:
                ax.text(left[i] + width / 2, i, f'{width:.1f}%', ha='center', va='center', fontsize=10)
        left = [a + b for a, b in zip(left, widths)]
    ax.set_yticks([0, 1], [f"Fungal entries (501)\n{grouped['ingroup']['representatives']:,} proteins",
                            f"Outgroups (25)\n{grouped['outgroup']['representatives']:,} proteins"])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel('Percentage of selected representative proteins')
    ax.set_title('Coding evidence and model availability across the full study')
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(loc='upper center', bbox_to_anchor=(0.43, -0.22), ncol=2, frameon=False, fontsize=9)
    fig.text(0.02, 0.02, 'Descriptive, protein-weighted source counts; no selection eligibility or predictor-confidence qualification.', fontsize=8)
    fig.subplots_adjust(left=0.24, right=0.97, top=0.86, bottom=0.34)
    for output in [args.png, args.pdf]:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=180)
    plt.close(fig)
    sources = {str(q): sha(q) for q in [args.completion, table, Path(__file__)]}
    outputs = {str(q): sha(q) for q in [args.png, args.pdf]}
    result = dict(status='complete_descriptive_full_coding_structure_source_coverage_figure',
                  checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526,
                  representative_proteins=5815847, grouped_counts=grouped,
                  source_hashes=sources, output_hashes=outputs, scientific_eligibility=False,
                  matplotlib_version=matplotlib.__version__, gpu=False, new_predictions=0,
                  scope='All526taxon rows and every source cross-count independently recounted for '
                        'a descriptive protein-weighted figure. No taxon-balanced estimate, biological '
                        'replication, phylogenetic effect, independent translation, selection or '
                        'predictor confidence/accuracy qualification.')
    with args.receipt.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
