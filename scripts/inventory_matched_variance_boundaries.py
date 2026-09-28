"""Census exact zero components in the original full fitted grid, without filtering."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ancestral_chain_attempt import sha, write_json


def main():
    # Original export predates the separate, still-running refinement overlay.
    root = Path('results/structural_comparisons/full-working-model-grid-export-20260927-v1')
    receipt_path = root / 'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'complete_full_working_model_grid_export'
    path = root / 'unique_fits.parquet'
    assert sha(path) == receipt['artifacts'][path.name]
    frame = pd.read_parquet(path)
    assert len(frame) == receipt['unique_fits'] == 144040
    assert not frame.duplicated(['fit_input_id', 'tree']).any()
    assert frame.fit_input_id.nunique() == 28808 and frame.tree.nunique() == 5
    names = ['background', 'family_component', 'species']
    ratios = frame[['variance_ratio_' + name for name in names]].to_numpy(float)
    finite = np.isfinite(ratios).all(axis=1)
    assert np.all(ratios[finite] >= 0)
    zero = ratios == 0
    patterns = [','.join(name for name, flag in zip(names, row) if flag) or 'none'
                for row in zero]
    patterns = np.where(finite, patterns, 'unavailable')
    selected = frame[['fit_input_id', 'tree', 'records', 'numerical_review_required',
                      'source_fit_sha256']].copy()
    selected['exact_zero_components'] = patterns
    selected['dense_float64_matrix_gib'] = selected.records.astype(float)**2 * 8 / 2**30
    counts = selected.groupby(['tree', 'numerical_review_required', 'exact_zero_components'],
                              dropna=False).size().reset_index(name='fits')
    assert counts.fits.sum() == len(frame)
    out = Path('results/model_validation/matched-variance-boundaries-20260928-v1')
    out.mkdir(parents=True, exist_ok=False)
    selected.to_parquet(out / 'fit_inventory.parquet', index=False)
    counts.to_csv(out / 'boundary_counts.tsv', sep='\t', index=False)
    pd.testing.assert_frame_equal(pd.read_parquet(out / 'fit_inventory.parquet'), selected)
    result = dict(status='full_original_fit_boundary_census', fits=len(frame),
        finite_variance_fits=int(finite.sum()), unavailable_variance_fits=int((~finite).sum()),
        any_exact_zero_fits=int((finite & zero.any(axis=1)).sum()),
        all_positive_fits=int((finite & ~zero.any(axis=1)).sum()),
        component_exact_zero_counts={name: int((finite & zero[:, i]).sum())
                                     for i, name in enumerate(names)},
        numerical_review_fits=int(frame.numerical_review_required.sum()),
        largest_single_dense_matrix_gib=float(selected.dense_float64_matrix_gib.max()),
        source_receipt_sha256=sha(receipt_path), source_export_sha256=sha(path),
        script_sha256=sha(__file__), artifacts={p.name: sha(p) for p in out.iterdir()},
        scope='All original unique tree fits, including numerical-review cases. Exact zeros '
              'only, no near-zero threshold. Refinement overlay is not included. Boundary '
              'estimates are not automatically invalid; this census does not determine '
              'identifiability, information rank, interval coverage or biological effects. '
              'Dense memory is one float64 n-by-n matrix, not peak memory or runtime.')
    write_json(out / 'receipt.json', result)
    write_json(Path('metadata/matched_variance_boundaries_20260928.json'), result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
