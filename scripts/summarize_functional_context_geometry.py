"""Summarize all four coordinate definitions and state/partner strata."""
import csv
import hashlib
import json
import math
from pathlib import Path
import pandas as pd


def main():
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    root = Path('results/functional_sites/prediction-context-geometry-20260927-v1')
    receipt = json.loads((root / 'receipt.json').read_text())
    source = root / 'functional_context_geometry.tsv'
    assert receipt['status'] == 'complete_functional_context_coordinate_comparison'
    assert sha(source) == receipt['artifacts'][source.name]
    data = pd.read_csv(source, sep='\t')
    keys = ['context_definition', 'paired_state_comparison', 'partner_comparison']
    summary = data.groupby(keys).agg(
        residues=('ca_rmsd_angstrom', 'size'), rmsd_min=('ca_rmsd_angstrom', 'min'),
        rmsd_median=('ca_rmsd_angstrom', 'median'), rmsd_max=('ca_rmsd_angstrom', 'max'),
        median_pair_distance_absolute_difference=('pair_distance_mean_absolute_difference_angstrom', 'median'))
    target = Path('metadata/functional_context_geometry_summary_20260927.tsv')
    assert not target.exists()
    summary.reset_index().to_csv(target, sep='\t', index=False)
    original = list(csv.DictReader(source.open(), delimiter='\t'))
    for row in csv.DictReader(target.open(), delimiter='\t'):
        chosen = [r for r in original if all(r[k] == row[k] for k in keys)]
        values = sorted(float(r['ca_rmsd_angstrom']) for r in chosen)
        distances = sorted(float(r['pair_distance_mean_absolute_difference_angstrom']) for r in chosen)
        n = len(values)
        assert int(row['residues']) == n
        median = lambda x: (x[(n-1)//2] + x[n//2])/2
        expected = dict(rmsd_min=values[0], rmsd_median=median(values), rmsd_max=values[-1],
                        median_pair_distance_absolute_difference=median(distances))
        assert all(math.isclose(float(row[k]), value, rel_tol=1e-12, abs_tol=1e-14) for k, value in expected.items())
    assert len(summary) == 16 and summary.residues.sum() == 600
    result = dict(status='passed_complete_functional_context_geometry_summary_readback',
        groups=16, comparisons=600, source_receipt_sha256=sha(root / 'receipt.json'),
        source_table_sha256=sha(source), output_sha256=sha(target), script_sha256=sha(__file__),
        scope='Every context definition and state/partner stratum retained. All minima, medians, maxima and counts checked by scalar sorting; descriptive dependent observations, not uncertainty intervals or inferential tests.')
    with Path('metadata/functional_context_geometry_summary_readback_20260927.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    with Path('metadata/functional_context_geometry_receipt_20260927.json').open('x') as handle:
        handle.write((root / 'receipt.json').read_text())
    print('Checked all 16 strata and 600 comparisons.')


if __name__ == '__main__':
    main()
