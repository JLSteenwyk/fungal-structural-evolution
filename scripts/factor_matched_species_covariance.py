#!/usr/bin/env python3
"""Factor complete working species covariance on its estimable design space."""
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from screen_duplication_domain_alignment_coverage import sha


def main():
    design_root = Path('results/phylogeny/matched-species-contrasts-20260927-v1')
    kernel_root = Path('results/phylogeny/species-distance-kernels-20260927-v1')
    out = Path('results/phylogeny/matched-species-covariance-factors-20260927-v1')
    pins = {}
    for root in [design_root, kernel_root]:
        pins[str(root / 'receipt.json')] = sha(root / 'receipt.json')
        for name,digest in json.loads((root / 'receipt.json').read_text())['artifacts'].items():
            assert sha(root / name) == digest
            pins[str(root / name)] = digest
    audit_path = Path('metadata/matched_species_contrast_readback_20260927.json')
    audit = json.loads(audit_path.read_text())
    assert audit['status'] == 'passed_full_matched_species_contrast_readback'
    assert audit['source_receipt_sha256'] == sha(design_root / 'receipt.json')
    pins[str(audit_path)] = sha(audit_path)
    design = sparse.load_npz(design_root / 'species_contrast_design.npz')
    u,s,vt = np.linalg.svd(design.toarray(), full_matrices=False)
    tolerance = max(design.shape) * np.finfo(float).eps * s[0]
    rank = int((s > tolerance).sum())
    assert rank == audit['contrast_rank'] == 242
    q = u[:, :rank]
    reduced = s[:rank, None] * vt[:rank]
    np.testing.assert_allclose(q @ reduced, design.toarray(), rtol=1e-10, atol=1e-12)
    out.mkdir(exist_ok=False)
    np.save(out / 'pattern_basis.npy', q)
    np.save(out / 'species_projection.npy', reduced)
    (out / 'patterns.tsv').write_bytes((design_root / 'patterns.tsv').read_bytes())
    (out / 'taxa.json').write_bytes((design_root / 'taxa.json').read_bytes())
    summaries = {}
    for path in sorted(kernel_root.glob('*.npz')):
        with np.load(path) as values:
            kernel = values['centered_kernel']
        core = reduced @ kernel @ reduced.T
        assert np.max(abs(core-core.T)) < 1e-10
        core = (core+core.T)/2
        # Failure here is a model/numerical diagnostic, never an invitation to add jitter.
        lower = np.linalg.cholesky(core)
        factor = q @ lower
        np.savez_compressed(out / (path.stem + '.npz'), factor=factor, reduced_covariance=core, cholesky=lower)
        maximum = 0.
        checked = 0
        for start in range(0, design.shape[0], 128):
            end = min(start+128, design.shape[0])
            exact = (design[start:end] @ kernel) @ design.T
            recovered = factor[start:end] @ factor.T
            np.testing.assert_allclose(recovered, exact, atol=1e-11, rtol=1e-10)
            maximum = max(maximum, float(np.max(abs(recovered-exact))))
            checked += exact.size
        summaries[path.stem] = dict(covariance_entries_checked=checked, maximum_absolute_error=maximum,
                                   minimum_core_eigenvalue=float(np.linalg.eigvalsh(core)[0]))
        print('Factored and checked all covariance entries:', path.stem, flush=True)
    for path,digest in pins.items():
        assert sha(path) == digest, path
    receipt = dict(status='complete_matched_species_covariance_factors_pending_tree_readback',
                   patterns=design.shape[0], species_columns=design.shape[1], rank=rank,
                   svd_rank_tolerance=float(tolerance), trees=summaries,
                   pins=pins, artifacts={p.name:sha(p) for p in out.iterdir()}, script_sha256=sha(__file__),
                   scope='Working W C W-transpose represented as F F-transpose on the measured design row space. All pattern-pair covariance entries checked. No full species inverse, eigenvalue clipping, diagonal jitter, fitted variance, or biological inference. Same pre-qualification design and five tree alternatives; per-setting row selection remains necessary.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['pins','artifacts']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
