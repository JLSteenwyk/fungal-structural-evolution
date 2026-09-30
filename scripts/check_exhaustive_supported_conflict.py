#!/usr/bin/env python3
"""Exercise independent full conflict maxima and rejection of a false stored maximum."""
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from screen_duplication_alignment_reuse import sha


def main():
    with tempfile.TemporaryDirectory() as name:
        root = Path(name)
        snapshot, audit, trees = [root / x for x in ['snapshot', 'audit', 'trees']]
        for path in [snapshot, audit, trees]:
            path.mkdir()
        guides = [root / 'guide1', root / 'guide2']
        full = '((A:1,B:1):1,(C:1,D:1):1,E:1);'
        alternative = '((A:1,C:1)90:1,(B:1,D:1)95:1,E:1);'
        for guide, tree in zip(guides, [full, alternative]):
            guide.mkdir()
            (guide / 'guide.treefile').write_text(tree)
            (guide / 'receipt.json').write_text(json.dumps(dict(
                status='passed_full_species_guide_readback',
                artifacts={'guide.treefile': sha(guide / 'guide.treefile')})))
        definitions = [('m1', alternative, [('A;C', '90'), ('B;D', '95')]),
                       ('m2', '((A:1,B:1)50:1,C:1,E:1);', [('A;B', '50')])]
        inputs, support = [], []
        for marker, tree, splits in definitions:
            folder = trees / marker
            folder.mkdir()
            path = folder / 'tree.treefile'
            path.write_text(tree)
            inputs.append(dict(marker=marker, tree_sha256=sha(path)))
            support.extend(dict(marker=marker, smaller_side_taxa=side, sh_alrt_percent=score)
                           for side, score in splits)
        with (snapshot / 'branch_support.tsv').open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(support[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(support)
        (snapshot / 'receipt.json').write_text(json.dumps(dict(
            completed_markers=2, planned_markers=2, pending_markers=[], inputs=inputs,
            artifacts={'branch_support.tsv': sha(snapshot / 'branch_support.tsv')})))
        (audit / 'receipt.json').write_text(json.dumps(dict(
            status='passed_all_snapshot_input_and_graph_split_readbacks',
            snapshot_receipt_sha256=sha(snapshot / 'receipt.json'), artifacts={})))
        diagnostic = root / 'diagnostic'
        subprocess.run([sys.executable, 'scripts/assess_supported_marker_guide_conflict.py',
                        '--snapshot', str(snapshot), '--readback', str(audit), '--trees', str(trees),
                        '--guides', *map(str, guides), '--output', str(diagnostic)],
                       check=True, capture_output=True, text=True)
        command = [sys.executable, 'scripts/readback_supported_marker_conflict_v2.py',
                   '--snapshot', str(snapshot), '--trees', str(trees), '--diagnostic', str(diagnostic)]
        proof = root / 'proof.json'
        subprocess.run(command + ['--output', str(proof)], check=True, capture_output=True, text=True)
        record = json.loads(proof.read_text())
        assert record['grid_status_and_witness_rows_checked'] == 16
        assert record['independent_full_set_conflict_maxima_checked'] == 8
        table = diagnostic / 'marker_guide_conflict.tsv'
        with table.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        assert {r['status'] for r in rows} == {
            'supported_concordance', 'supported_conflict', 'unresolved_at_support_cutoff',
            'uninformative_taxon_coverage'}
        original = [dict(row) for row in rows]
        for row in rows:
            if (row['marker'], row['guide'], row['full_split_taxa'], row['sh_alrt_cutoff']) == ('m1', 'guide1', 'A;B', '95'):
                row['maximum_conflicting_sh_alrt'] = '90'
                row['maximum_conflict_witness_split'] = 'A;C'
                row['status'] = 'unresolved_at_support_cutoff'
        with table.open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        receipt_path = diagnostic / 'receipt.json'
        receipt = json.loads(receipt_path.read_text())
        receipt['artifacts'][table.name] = sha(table)
        receipt_path.write_text(json.dumps(receipt))
        failure = subprocess.run(command + ['--output', str(root / 'bad-cutoff-proof.json')],
                                 capture_output=True, text=True)
        assert failure.returncode != 0 and 'nunique(dropna=False)' in failure.stderr
        assert not (root / 'bad-cutoff-proof.json').exists()
        rows = original
        # Rehash the edited artifact so rejection must arise from reconstruction.
        for row in rows:
            if row['marker'] == 'm1' and row['guide'] == 'guide1' and row['full_split_taxa'] == 'A;B':
                row['maximum_conflicting_sh_alrt'] = '90'
                row['maximum_conflict_witness_split'] = 'A;C'
                if row['sh_alrt_cutoff'] == '95':
                    row['status'] = 'unresolved_at_support_cutoff'
        with table.open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        receipt_path = diagnostic / 'receipt.json'
        receipt = json.loads(receipt_path.read_text())
        receipt['artifacts'][table.name] = sha(table)
        receipt_path.write_text(json.dumps(receipt))
        failure = subprocess.run(command + ['--output', str(root / 'bad-proof.json')],
                                 capture_output=True, text=True)
        assert failure.returncode != 0 and 'assert (max(scores)' in failure.stderr
        assert not (root / 'bad-proof.json').exists()
    print('All four support/coverage states, exhaustive maxima, complete grid, cutoff-consistency and rehashed false-maximum rejection passed.')


if __name__ == '__main__':
    main()
